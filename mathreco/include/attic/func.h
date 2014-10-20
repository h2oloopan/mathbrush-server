#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_


namespace scg {


template <typename R, typename T>
class bound_mem_fun_t {
public:
	typedef R ret_type;

private:
	typedef R (T::*FnType)();

	FnType fn;
	T *obj;

public:
	bound_mem_fun_t(FnType fn_, T *obj_) : fn(fn_), obj(obj_) { }
	R operator()() {
		return (obj->*fn)();
	}
};

template <typename R, typename T>
class const_bound_mem_fun_t {
public:
	typedef R ret_type;

private:
	typedef R (T::*FnType)() const;

	FnType fn;
	const T *obj;

public:
	const_bound_mem_fun_t(FnType fn_, const T *obj_) : fn(fn_), obj(obj_) { }
	R operator()() const {
		return (obj->*fn)();
	}
};

template<typename R, typename T>
bound_mem_fun_t<R, T>
bound_mem_fun(R (T::*fn)(), T *obj) {
	return bound_mem_fun_t<R, T>(fn, obj);
}
template<typename R, typename T>
const_bound_mem_fun_t<R, T>
const_bound_mem_fun(R (T::*fn)() const, const T *obj) {
	return const_bound_mem_fun_t<R, T>(fn, obj);
}


template <typename R, typename A, typename T>
class bound_mem_fun1_t {
public:
	typedef R ret_type;

private:
	typedef R (T::*FnType)(A);

	FnType fn;
	T *obj;

public:
	bound_mem_fun1_t(FnType fn_, T *obj_) : fn(fn_), obj(obj_) { }
	R operator()(A a) const {
		return (obj->*fn)(a);
	}
};

template <typename R, typename A, typename T>
class const_bound_mem_fun1_t {
public:
	typedef R ret_type;

private:
	typedef R (T::*FnType)(A) const;

	FnType fn;
	const T *obj;

public:
	const_bound_mem_fun1_t(FnType fn_, const T *obj_) : fn(fn_), obj(obj_) { }
	R operator()(A a) const {
		return (obj->*fn)(a);
	}
};

template<typename R, typename A, typename T>
bound_mem_fun1_t<R, A, T>
bound_mem_fun1(R (T::*fn)(A), T *obj) {
	return bound_mem_fun1_t<R, A, T>(fn, obj);
}
template<typename R, typename A, typename T>
const_bound_mem_fun1_t<R, A, T>
const_bound_mem_fun1(R (T::*fn)(A) const, const T *obj) {
	return const_bound_mem_fun1_t<R, A, T>(fn, obj);
}

}


#endif
