#include "analysis7b.c"

void job_wrapper( const char* input_file_list="", int job_index = 0 ) {
   printf("\n\n Inside job_wrapper:  %s, %d\n\n", input_file_list, job_index ) ;
   analysis7b a( input_file_list, job_index ) ;
   a.Loop(0) ;
   printf("\n\n Done in job_wrapper.\n\n") ;
}



