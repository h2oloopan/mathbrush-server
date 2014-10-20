#include "proc.h"
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdlib.h>

#define _WIN32_WINNT 0x0501
#define WIN32_LEAN_AND_MEAN
#include <windows.h>


namespace scg
{


struct proc_t
{
    LPVOID fiber;
};


proc_t *
proc_create_base()
{
	proc_t *p = new proc_t;
	if (!p) {
		return 0;
	}

	p->fiber = ConvertThreadToFiber(0);
	return p;
}


proc_t *
proc_create(proc_fn cb, void *data, unsigned stack_size)
{
	proc_t *p;
	
	assert(cb);
	
	p = new proc_t;

    p->fiber = CreateFiberEx(stack_size, stack_size, 0, (LPFIBER_START_ROUTINE)cb, data);
    
    if (!p->fiber) {
        delete p;
        return 0;
    }
    
	return p;
}

void
proc_destroy(proc_t *p)
{
	assert(p);
    DeleteFiber(p->fiber);
	delete p;
}

void
proc_switch(proc_t *p)
{
    SwitchToFiber(p->fiber);
}

}