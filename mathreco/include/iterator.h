#ifndef ITERATOR_H_
#define ITERATOR_H_


#include <cassert>


namespace scg
{


template <typename ContainerType, typename ContainerIteratorType, typename IteratorType>
class MultiContainerIterator
{
public:
    typedef typename IteratorType::pointer pointer;
    typedef typename IteratorType::reference reference;
    
    
public:
    MultiContainerIterator() : base_container(0) {}
    
    explicit MultiContainerIterator(ContainerType &cont)
      : base_container(&cont)
    {
        if (base_container->empty()) {
            base_container = 0;
        }
        else {
            base_iter = base_container->begin();
            iter = base_iter->begin();
        }
    }

    MultiContainerIterator(ContainerType &cont,
                           ContainerIteratorType base_start)
      : base_container(&cont), base_iter(base_start)
	 {
	 	if (base_iter != base_container->end()) {
			iter = base_iter->begin();
		}
	 }

    MultiContainerIterator(ContainerType &cont,
                           ContainerIteratorType base_start,
                           IteratorType iter_start)
      : base_container(&cont), base_iter(base_start), iter(iter_start) {}

    pointer operator->()
      { return iter.operator->(); }
    
    reference operator*()
      { return *iter; }

    const pointer operator->() const
      { return iter.operator->(); }
    
    reference operator*() const
      { return *iter; }
    
    MultiContainerIterator &operator++() // preincrement
    {
        assert(base_container);
        if (base_iter < base_container->end()) {
            ++iter;
            if (iter >= base_iter->end()) {
                ++base_iter;
                if (base_iter == base_container->end()) {
                    iter = base_container->begin()->end();
                }
                else {
                    iter = base_iter->begin();
                }
            }
        }
        return *this;
    }
    
    MultiContainerIterator operator++(int) // postincrement
    {
        MultiContainerIterator tmp(*this);
        ++(*this);
        return tmp;
    }
    
    MultiContainerIterator &operator--() // predecrement
    {
        assert(base_container);
        if (base_iter >= base_container->begin()) {
            --iter;
            if (iter <= base_iter->begin()) {
                --base_iter;
                iter = base_iter->end() - 1;
            }
        }
        return *this;
    }

    MultiContainerIterator operator--(int) // postdecrement
    {
        MultiContainerIterator tmp(*this);
        --(*this);
        return tmp;    
    }

    MultiContainerIterator &operator+=(int n)
    {   // TODO: improve this implementation
        while (n--) {
            ++(*this);
        }
        return *this;
    }
    
    MultiContainerIterator &operator-=(int n)
    {   // TODO: improve this implementation
        while (n--) {
            --(*this);
        }
        return *this;
    }
    
    MultiContainerIterator operator+(int rhs) const
    {
        MultiContainerIterator it(*this);
        it += rhs;
        return it;
    }
    
    MultiContainerIterator operator-(int rhs) const
    {
        MultiContainerIterator it(*this);
        it -= rhs;
        return it;
    }
    
    bool operator==(const MultiContainerIterator &rhs) const
    {
        return (!base_container && !rhs.base_container)
            || ((base_iter == rhs.base_iter) && ((base_iter == base_container->end()) || (iter == rhs.iter)));
    }

    bool operator!=(const MultiContainerIterator &rhs) const
    {
        if (!base_container && !rhs.base_container) {
            return false;
        }
        else if (base_container != rhs.base_container) {
            return true;
        }
        
        return (base_iter != rhs.base_iter) || ((base_iter != base_container->end()) && (iter != rhs.iter));
    }

    bool operator<=(const MultiContainerIterator &rhs) const
    {
        return (base_iter < rhs.base_iter) || ((base_iter == rhs.base_iter) && (iter <= rhs.base_iter));
    }

    bool operator>=(const MultiContainerIterator &rhs) const
    {
        return (base_iter > rhs.base_iter) || ((base_iter == rhs.base_iter) && (iter >= rhs.base_iter));
    }

    bool operator<(const MultiContainerIterator &rhs) const
    {
        return (base_iter < rhs.base_iter) || ((base_iter == rhs.base_iter) && (iter < rhs.base_iter));
    }

    bool operator>(const MultiContainerIterator &rhs) const
    {
        return (base_iter > rhs.base_iter) || ((base_iter == rhs.base_iter) && (iter > rhs.base_iter));
    }


private:
    ContainerType *base_container;
    ContainerIteratorType base_iter;
    IteratorType iter;
};



template <typename IteratorType, typename PatternType>
class SelectiveIterator
{
public:
    typedef typename IteratorType::pointer pointer;
    typedef typename IteratorType::reference reference;
    
public:
    SelectiveIterator() {}
    SelectiveIterator(IteratorType s, IteratorType e, const PatternType &pat) : start(s), end(e), iter(start), pattern(pat) {}
    
    pointer operator->()
    {
        return iter.operator->();
    }
    
    reference operator*()
    {
        return *iter;
    }
    
    SelectiveIterator &operator++() // preincrement
    {
        // use !(A == B) so that only operator== must be defined by clients
        while (!(*(++iter) == pattern) && !(iter == end));
        return *this;
    }
    
    SelectiveIterator operator++(int) // postincrement
    {
        SelectiveIterator tmp(*this);
        ++(*this);
        return tmp;
    }
    
    SelectiveIterator &operator--() // predecrement
    {
        while (*(--iter) != pattern && iter >= start);
        return *this;
    }

    SelectiveIterator operator--(int) // postdecrement
    {
        SelectiveIterator tmp(*this);
        --(*this);
        return tmp;    
    }

    bool operator==(const SelectiveIterator &rhs) const
    {
        return iter == rhs.iter;
    }

    bool operator!=(const SelectiveIterator &rhs) const
    {
        return iter != rhs.iter;
    }

    bool operator<=(const SelectiveIterator &rhs) const
    {
        return iter <= rhs.iter;
    }

    bool operator>=(const SelectiveIterator &rhs) const
    {
        return iter <= rhs.iter;
    }

    bool operator<(const SelectiveIterator &rhs) const
    {
        return iter < rhs.iter;
    }

    bool operator>(const SelectiveIterator &rhs) const
    {
        return iter > rhs.iter;
    }

    bool operator==(const IteratorType &rhs) const
    {
        return iter == rhs;
    }

    bool operator!=(const IteratorType &rhs) const
    {
        return iter != rhs;
    }

    bool operator<=(const IteratorType &rhs) const
    {
        return iter <= rhs;
    }

    bool operator>=(const IteratorType &rhs) const
    {
        return iter <= rhs;
    }

    bool operator<(const IteratorType &rhs) const
    {
        return iter < rhs;
    }

    bool operator>(const IteratorType &rhs) const
    {
        return iter > rhs;
    }

public:
    IteratorType start;
    IteratorType end;
    IteratorType iter;
    
    PatternType pattern;
};


}


#endif

