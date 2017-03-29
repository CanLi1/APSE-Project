#include "wrapper.h"
#include <stdio.h>
#include <assert.h>

FUNCGEN_DUMMY(V_vap_dummy);

int main(void) 
{
   arglist* al;
   double y = 0.0;
   y = V_vap_dummy(al);
   assert(y==-1000.0);
   return 0;
}