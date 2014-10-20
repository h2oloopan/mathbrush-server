#include "proc.h"
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdlib.h>

#include <sys/mman.h>
#include <sys/select.h>

#ifndef PAGESIZE
#define PAGESIZE 4096
#endif



namespace scg
{


struct proc_t
{
	proc_fn cb;
	void *dat;
	jmp_buf jmp;
	unsigned char *stack;
	int running;
};

static proc_t *curr_proc = NULL;

static void
proc_startup(int unused)
{
	if (!setjmp(curr_proc->jmp)) {
		curr_proc->running = 1;
		return;
	}

	curr_proc->cb(curr_proc->dat);
}

static int
proc_bootstrap(proc_t *p, unsigned stack_size)
{
	const static unsigned PROC_STACK_SIZE = PAGESIZE * 32;

	if (stack_size == 0) {
		stack_size = PROC_STACK_SIZE;
	}

	struct sigaltstack stack;
	struct sigaction sa, old_sa;
	sigset_t sigs, old_sigs;
	proc_t *curr = curr_proc;

	stack.ss_flags = 0;
	stack.ss_size = PROC_STACK_SIZE;// - sizeof(int);
	if (!(p->stack = new unsigned char[stack_size + PAGESIZE - 1])) {
		return -1;
	}
	stack.ss_sp = (void *)(((unsigned int)p->stack + PAGESIZE - 1) & ~(PAGESIZE - 1));
	sigaltstack(&stack, NULL);

	sigemptyset(&sigs);
	sigaddset(&sigs, SIGUSR1);
	sigprocmask(SIG_BLOCK, &sigs, &old_sigs);

	sa.sa_handler = &proc_startup;
	sa.sa_flags = SA_ONSTACK;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGUSR1, &sa, &old_sa);

	curr_proc = p;

	raise(SIGUSR1);

	sigfillset(&sigs);
	sigdelset(&sigs, SIGUSR1);
	while (!p->running) {
		sigsuspend(&sigs);
	}

	curr_proc = curr;

	sigaltstack(NULL, &stack);
	stack.ss_flags = SS_DISABLE;
	sigaltstack(&stack, NULL);
	sigaltstack(NULL, &stack);
	
	sigaction(SIGUSR1, &old_sa, NULL);
	sigprocmask(SIG_SETMASK, &old_sigs, NULL);

	return 0;
}


proc_t *
proc_create_base()
{
	proc_t *p = (proc_t *)malloc(sizeof(proc_t));
	if (!p) {
		return 0;
	}

	p->stack = 0;
	p->cb = 0;
	p->dat = 0;
	p->running = 1;
	if (!setjmp(p->jmp)) {
		proc_switch(p);
	}
	return p;
}


proc_t *
proc_create(proc_fn cb, void *data, unsigned stack_size)
{
	proc_t *p;
	
	assert(cb);
	
	p = new proc_t;

	p->cb = cb;
	p->dat = data;
	p->running = 0;
	
	if (proc_bootstrap(p, stack_size) < 0) {
		delete p;
		return NULL;
	}

	return p;
}

void
proc_destroy(proc_t *p)
{
	assert(p);
	delete[] p->stack;
	delete p;
}

void
proc_switch(proc_t *p)
{
	if (!curr_proc) {
		curr_proc = p;
		longjmp(p->jmp, 1);
	}
	else if (!setjmp(curr_proc->jmp)) {
		curr_proc = p;
		longjmp(p->jmp, 1);
	}
}


}
