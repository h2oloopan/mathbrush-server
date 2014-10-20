/*
 * functional.h
 *
 * This file defines extensions to STL's functional programming API.
 * In particular, it defines generic functors for the composition of
 * two functors.
 *
 */
 
#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_


#include <functional>


namespace scg
{


template <typename T>
void destroy(T *t) { delete t; }


template <typename OuterFn, typename InnerFn>
struct unary_binary_composition_t : public std::binary_function<typename InnerFn::first_argument_type,
                                                                typename InnerFn::second_argument_type,
                                                                typename OuterFn::result_type>
{
public:
    unary_binary_composition_t(OuterFn F, InnerFn G) : f(F), g(G) {}
    
    typename OuterFn::result_type operator()(typename InnerFn::first_argument_type a,
                                             typename InnerFn::second_argument_type b) const
    {
        return f(g(a, b));
    }
    
    
private:
    OuterFn f;
    InnerFn g;
};


template <typename OuterFn, typename InnerFn>
unary_binary_composition_t<OuterFn, InnerFn>
unary_binary_composition(OuterFn f, InnerFn g)
{
    return unary_binary_composition_t<OuterFn, InnerFn>(f, g);
}


template <typename OuterFn, typename InnerFn>
struct unary_unary_composition_t : public std::unary_function<typename InnerFn::argument_type, typename OuterFn::result_type>
{
public:
    unary_unary_composition_t(OuterFn F, InnerFn G) : f(F), g(G) {}
    
    typename OuterFn::result_type operator()(typename InnerFn::argument_type a) const
    {
        return f(g(a));
    }

private:
    OuterFn f;
    InnerFn g;
};


template <typename OuterFn, typename InnerFn>
unary_unary_composition_t<OuterFn, InnerFn>
unary_unary_composition(OuterFn f, InnerFn g)
{
    return unary_unary_composition_t<OuterFn, InnerFn>(f, g);
}


template <typename OuterFn, typename InnerFn>
struct binary_unary1_composition_t : public std::binary_function<typename OuterFn::first_argument_type,
                                                                 typename InnerFn::argument_type,
                                                                 typename OuterFn::result_type>
{
    binary_unary1_composition_t(OuterFn F, InnerFn G) : f(F), g(G) {}
    
    typename OuterFn::result_type operator()(typename OuterFn::first_argument_type a1,
                                             typename InnerFn::argument_type a2) const
    {
        return f(g(a1), a2);
    }

private:
    OuterFn f;
    InnerFn g;
};


template <typename OuterFn, typename InnerFn>
binary_unary1_composition_t<OuterFn, InnerFn>
binary_unary1_composition(OuterFn f, InnerFn g)
{
    return binary_unary1_composition_t<OuterFn, InnerFn>(f, g);
}


template <typename OuterFn, typename InnerFn>
struct binary_unary2_composition_t : public std::binary_function<typename OuterFn::first_argument_type,
                                                                 typename InnerFn::argument_type,
                                                                 typename OuterFn::result_type>
{
    binary_unary2_composition_t(OuterFn F, InnerFn G) : f(F), g(G) {}
    
    typename OuterFn::result_type operator()(typename OuterFn::first_argument_type a1,
                                             typename InnerFn::argument_type a2) const
    {
        return f(a1, g(a2));
    }

private:
    OuterFn f;
    InnerFn g;
};


template <typename OuterFn, typename InnerFn>
binary_unary2_composition_t<OuterFn, InnerFn>
binary_unary2_composition(OuterFn f, InnerFn g)
{
    return binary_unary2_composition_t<OuterFn, InnerFn>(f, g);
}


}


#endif

