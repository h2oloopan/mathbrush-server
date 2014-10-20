#ifndef PROC_H_
#define PROC_H_


namespace scg
{


struct proc_t;

typedef void (*proc_fn)(void *);


proc_t * proc_create_base();
proc_t * proc_create(proc_fn cb, void *data, unsigned stack_size);
void proc_destroy(proc_t *p);
void proc_switch(proc_t *p);


}


#endif
