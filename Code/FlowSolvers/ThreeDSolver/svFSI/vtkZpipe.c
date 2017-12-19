
   #include <stdio.h>
   #include <string.h>
   #include <assert.h>
   #include "simvascular_zlib.h"

#ifdef SV_WRAP_FORTRAN_IN_CAPS_NO_UNDERSCORE
  #define defzlibdata_ DEFZLIBDATA
#elif SV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  #define defzlibdata_ defzlibdata_
#endif

#ifdef SV_WRAP_FORTRAN_IN_CAPS_NO_UNDERSCORE
  #define infzlibdata_ INFZLIBDATA
#elif SV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  #define infzlibdata_ infzlibdata_
#endif

   int inf ( unsigned char* in, int nin, unsigned char* out, int nout );
   int def ( unsigned char* in, int nin, unsigned char* out, int* nout, int level);	
   
   void zerr(int ret)
   {
      fputs("zpipe: ", stderr);
      switch (ret) {
      case Z_ERRNO:
         if (ferror(stdin))
            fputs("Error reading stdin\n", stderr);
         if (ferror(stdout))
            fputs("Error writing stdout\n", stderr);
         break;
      case Z_STREAM_ERROR:
         fputs("Invalid compression level\n", stderr);
         break;
      case Z_DATA_ERROR:
         fputs("Invalid or incomplete deflate data\n", stderr);
         break;
      case Z_MEM_ERROR:
         fputs("Out of memory\n", stderr);
         break;
      case Z_VERSION_ERROR:
         fputs("zlib version mismatch!\n", stderr);
      }
   }
      
   extern int defzlibdata_( unsigned char* in, int* m, unsigned char* out, int* n, int* plev, int* ierr ) {
      int ret;
      int nin = *m;
      int level = *plev;
      
      *ierr = 0;
      ret = def( in, nin, out, n, level);
      if ( ret!=Z_OK ) {
         zerr(ret);
         *ierr = -1;
      }
      return ret;
   }
   
   extern int infzlibdata_( unsigned char* in, int* m, unsigned char* out, int* n, int* ierr ) {
      int ret;
      int nin = *m;
      int nout = *n;
      
      *ierr = 0;
      ret = inf( in, nin, out, nout);
      if ( ret!=Z_OK ) {
         zerr(ret);
         *ierr = -1;
      }
      return ret;
   }

   int def ( unsigned char* in, int nin, unsigned char* out, int* nout, int level ) {
      int ret, flush;
      char c;
      
      z_stream strm;
      strm.zalloc = 0;
      strm.zfree = 0;
      ret = deflateInit(&strm, level);
      if (ret != Z_OK)
         return ret;

      strm.avail_in = nin;
      strm.next_in = in;
      strm.avail_out = *nout;
      strm.next_out = out;
      
      while (strm.avail_in != 0)
      {
          ret = deflate(&strm, flush);
          assert(ret != Z_STREAM_ERROR);
          if ( strm.avail_out == 0 )
          {
              strm.avail_out = *nout;
              strm.next_out = out;
          }
      }
      
      ret = Z_OK;
      while(ret == Z_OK)
      {
          flush = Z_FINISH;
          if ( strm.avail_out == 0 )
          {
              strm.avail_out = *nout;
              strm.next_out = out;
          }
          ret = deflate(&strm, flush);
	  }
      assert(ret == Z_STREAM_END);
      
      *nout = *nout - strm.avail_out;
      (void)deflateEnd(&strm);
      return Z_OK;
      
   }
   
   int inf ( unsigned char* in, int nin, unsigned char* out, int nout ) {
      int ret;
      unsigned have;
      z_stream strm;
      
       /* allocate inflate state */
      strm.zalloc = Z_NULL;
      strm.zfree = Z_NULL;
      strm.opaque = Z_NULL;
      strm.avail_in = 0;
      strm.next_in = Z_NULL;
      ret = inflateInit(&strm);
      if (ret != Z_OK)
         return ret;

      strm.avail_in = nin;
      strm.next_in = in;

      /* run inflate() on input until output buffer not full */
      do {
         strm.avail_out = nout;
         strm.next_out = out;
         ret = inflate(&strm, Z_NO_FLUSH);
         assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
         switch (ret) {
         case Z_NEED_DICT:
            fprintf(stdout,"zpipe: need dictionary..\n");
            ret = Z_DATA_ERROR;     /* and fall through */
         case Z_DATA_ERROR:
            fprintf(stdout,"zpipe: data error..\n");
         case Z_MEM_ERROR:
            (void)inflateEnd(&strm);
            return ret;
         }
         have = nout - strm.avail_out;
      } while (strm.avail_out == 0);

      /* clean up and return */
      (void)inflateEnd(&strm);
      return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
   }

