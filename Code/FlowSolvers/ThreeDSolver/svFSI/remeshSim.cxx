#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#ifndef __MACH__
    #include "MeshSim.h"
    #include "SimExport.h"
    #include "SimUtil.h"
    #include "SimDiscrete.h"
    #include "SimMessages.h"
    #include "SimError.h"
    #include "SimErrorCodes.h"
    #include "SimMeshingErrorCodes.h"
    #include "SimDiscreteErrorCodes.h"
    #include "SimExportErrorCodes.h"
    
    void messageHandler(int type, const char *msg);
    void progressHandler(const char *what, int level, int startVal,  
                     int endVal, int currentVal, void *);
    void printModelStats(pGModel model);
    void writeNodesAndElements(pMesh, FILE *);
#endif

extern"C" {
   void remesh3d_meshsim_(const int* nPoints, const int* nFacets, const double* pointList, \
                    const int* facetList, const double* params, int* pOK)
   {
      FILE * fptr;
      #ifdef __MACH__
	       printf("ERROR: Cannot use MeshSim library on macOS\n");
	       *pOK = -1;
	       return;
      #else
         try {
            *pOK = 0;
            char fName[250];
            const int nDim = 3;
            const int nodesPerElem = 4;
            double maxRadRatio;
            double minDihedAng;
            double maxEdgeSize;
            double* coords = new double[nDim*(*nPoints)];
            int* elementData = new int[nodesPerElem*(*nFacets)];
            int* elementType = new int[*nFacets];

            for (int i=0; i < *nPoints; i++) {
                for (int j=0; j < nDim; j++) {
                    coords[nDim*i+j] = *pointList;
                    ++pointList;
                }
            }
        
            for (int i = 0; i < *nFacets; i++) {
                for (int j=0; j < nodesPerElem; j++) {
                    elementData[nodesPerElem*i+j] = (*facetList)-1;
                    ++facetList;
                }
                elementType[i] = 5;
            }
        
            maxRadRatio = *params;
            minDihedAng = *(++params);
            maxEdgeSize = *(++params);

            strcpy(fName, SIM_INC);
            strcat(fName,"/meshsim-license.dat");
            Sim_readLicenseFile(fName);

            pMesh oldMesh;
            pVertex* vReturn = new pVertex[*nPoints];
            pEntity* eReturn = new pEntity[*nFacets];
        
            // Initialize Meshsim
            Sim_logOn("meshsim.log");
            MS_init();
            SimDiscrete_start(0);
            Sim_setMessageHandler(messageHandler);
            pProgress progress = Progress_new();
            Progress_setCallback(progress, progressHandler);
        
            oldMesh = M_new(0,0);
            if(M_importFromData(oldMesh,*nPoints,coords,*nFacets,elementType,elementData,vReturn,eReturn,progress))
            {
                printf("error\n") ;
                M_release(oldMesh);
                *pOK = -1;
                return;
            }
        
            delete[] coords;
            delete[] elementData;
            delete[] elementType;

            // Check the input mesh for intersections
            // Note: This call must occur before the Discrete Model is created
            if(MS_checkMeshIntersections(oldMesh,0,progress))
            {
                printf("  WARNING!  There are intersections in the input mesh\n");
                M_release(oldMesh);
                *pOK = -1;
                return;
            }

            // We will use the geometry such that DM_findEdgeByFaceNormals will not 
            // cause issue and we should see all the model vertices labeled as well.
            // First, mark the vertices
            for(int i = 0; i < *nPoints; i++) {
                V_markAsBoundary(vReturn[i],i);
            }

            // The first 10 faces are in group 100. They are 100, 101, 102, 103, 104
            // The second 10 faces are in group 200. They are 200, 201, 202, 203, 204
            int increment = 0;
            int step = 0;
            for(int i = 0; i < *nFacets; i++)
            {
                if(!(i%2)) step++;
                if (i == 0 || i == 10 || i == 20 || 
                    i == 30 || i == 40 || i == 50) {
                    increment++;
                    step = 0;
                }
                //int faceValue = 100*increment+step;
                F_markAsBoundary((pFace)eReturn[i], i); //faceValue);
            }

            // Create the Discrete Model
            pDiscreteModel model = DM_createFromMesh(oldMesh, 0, progress);

            // Check for model error
            if(!model)
            {
                printf("  Error creating Discrete model from mesh\n");
                M_release(oldMesh);
                *pOK = -1;
                return;
            }

            DM_findEdgesByFaceNormals(model, 0, progress); 
            DM_eliminateDanglingEdges(model, progress);
            if( DM_completeTopology(model, progress) ) {
                printf("  Error completing Discrete model topology\n");
                M_release(oldMesh);
                GM_release(model);
                *pOK = -1;
                return;
            }

            pModelItem modelDomain = GM_domain(model);
            pMesh newMesh = DM_getMeshCopy(model,0 );   
            pACase meshCase = MS_newMeshCase(model);

            MS_setMeshSize(meshCase,modelDomain,1,maxEdgeSize,0);
            //pVolumeMeshImprover VMI = VolumeMeshImprover_new (newMesh);
            //VolumeMeshImprover_setModifySurface(VMI,0);
            pVolumeMesher volumeMesher = VolumeMesher_new (meshCase, newMesh);
            VolumeMesher_setModifySurface(volumeMesher,0);
            VolumeMesher_setOptimization(volumeMesher,0);
            VolumeMesher_setSmoothing(volumeMesher,0);
            VolumeMesher_execute(volumeMesher, progress);

            SimDiscrete_stop(0);

            SimExport_start();
            pExporter exporter = Exporter_new(); // Instantiate Exporter

            // Load Nastran pattern file
            strcpy(fName, SIM_INC);
            strcat(fName,"/amesh.sxp");
            pCompiledPattern pattern = CompiledPattern_createFromFile(fName);
            Exporter_setInputs(exporter, newMesh, 0);
            Exporter_setOutputs(exporter, "new-vol-mesh.node", "");
            Exporter_executeCompiledPattern(exporter, pattern, progress);

            CompiledPattern_release(pattern);
            Exporter_delete(exporter);

            //Repeat for aconn
            exporter = Exporter_new();
            strcpy(fName, SIM_INC);
            strcat(fName,"/aconn.sxp");
            pattern = CompiledPattern_createFromFile(fName);
            Exporter_setInputs(exporter, newMesh, 0);
            Exporter_setOutputs(exporter, "new-vol-mesh.ele", "");
            Exporter_executeCompiledPattern(exporter, pattern, progress);
            printf("Export complete\n");

            //Repeat for abound connectivity
            exporter = Exporter_new();
            strcpy(fName, SIM_INC);
            strcat(fName,"/abound_conn.sxp");
            pattern = CompiledPattern_createFromFile(fName);
            Exporter_setInputs(exporter, newMesh, 0);
            Exporter_setOutputs(exporter, "new-surf-mesh.ele", "");
            Exporter_executeCompiledPattern(exporter, pattern, progress);
            printf("Export complete\n");

            // cleanup
            delete[] eReturn;
            delete[] vReturn;
            M_release(newMesh);
            GM_release(model);
            Progress_delete(progress);

            // shutdown
            SimExport_stop();
            MS_exit();
            Sim_logOff();
            Sim_unregisterAllKeys();

            return;

         }
         catch (pSimError err) {
            printf("Simmetrix error caught:\n");
            printf("  Error code: %d\n",SimError_code(err));
            printf("  Error string: %s\n",SimError_toString(err));
            SimError_delete(err);
            return;
         } 
         catch (...) {
            printf("Unhandled exception caught\n");
            return;
         }
      #endif
    }
}

#ifndef __MACH__
    void printModelStats(pGModel model) {

        printf("  Number of model vertices: %d\n",GM_numVertices(model));
        printf("    Vertex tags: ");
        GVIter modelVertices = GM_vertexIter(model);  // initialize the iterator
        pGVertex modelVertex;
        while ( ( modelVertex = GVIter_next(modelVertices) ) ) { // get next vertex
            printf("%d ",GEN_tag(modelVertex));
        }
        GVIter_delete(modelVertices);
        printf("\n");
        printf("  Number of model edges: %d\n",GM_numEdges(model));
        printf("  Number of model faces: %d\n",GM_numFaces(model));
        printf("    Face tags: ");
        GFIter modelFaces = GM_faceIter(model);  // initialize the iterator
        pGFace modelFace;
        while ( ( modelFace = GFIter_next(modelFaces) ) ) { // get next face
            printf("%d ",GEN_tag(modelFace));
        }
        GFIter_delete(modelFaces);
        printf("\n");
        printf("  Number of model regions: %d\n",GM_numRegions(model));
    }

    void messageHandler(int type, const char *msg) {

        switch (type) {
            case Sim_InfoMsg:
                printf("Info: %s\n",msg);
                break;
            case Sim_DebugMsg:
                printf("Debug: %s\n",msg);
                break;
            case Sim_WarningMsg:
                printf("Warning: %s\n",msg);
                break;
            case Sim_ErrorMsg:
                printf("Error: %s\n",msg);
                break;
        }
        return;
    }

    void progressHandler(const char *what, int level, int startVal, 
                     int endVal, int currentVal, void *) {
  
        printf("Progress: %s, level: %d, startVal: %d, endVal: %d, currentVal: %d\n",
         what,level,startVal,endVal,currentVal);
        return;
    }

    void writeNodesAndElements(pMesh mesh, FILE *file) {

        VIter vertices;
        FIter faces;
        pVertex vertex;
        pFace face;
        double xyz[3];
        pPList faceVerts;
        int i,j;

        vertices = M_vertexIter(mesh);
        i=0;
        while ( ( vertex = VIter_next(vertices) ) ) {
            
            EN_setID((pEntity)vertex,i++);
            V_coord(vertex,xyz);
            fprintf(file,"%.9f %.9f %.9f\n",xyz[0],xyz[1],xyz[2]);
        }
        VIter_delete(vertices);

        faces = M_faceIter(mesh); //cnt=0;
        while ( ( face = FIter_next(faces) ) ) {
            
            faceVerts = F_vertices(face, 1);
            for (j=0; j<3; j++) {
                vertex = (pVertex)PList_item(faceVerts,j);
                fprintf(file,"%d ",EN_id((pEntity)vertex)+1);
            }
            fprintf(file,"\n");
            PList_delete(faceVerts);
        }
        FIter_delete(faces);

    }
#endif
