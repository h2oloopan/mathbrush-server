#ifndef NODE_H_
#define NODE_H_


#include <iterator>
#include <list>
#include <string>
#include <vector>

#include "group.h"
#include "match.h"

namespace scg
{

class ExpressionNode;

struct ExpressionLink
{
    ExpressionLink(ExpressionNode *from_, ExpressionNode *to_, const std::string &name_, double parm1_)
        : from(from_), to(to_), name(name_), parm1(parm1_) {}
    
    ExpressionNode *from;
    ExpressionNode *to;
    std::string name;
    double parm1;
};


class LinkCollection
{
public:
    class iterator : public std::iterator_traits<ExpressionLink *>
    {
    public:
        iterator(std::vector<ExpressionLink *> *coll_, std::vector<ExpressionLink *>::iterator start) : coll(coll_), it(start) {}
        
        iterator &operator=(const iterator &i)
        {
            if (&i == this) {
                return *this;
            }
            
            coll = i.coll;
            it = i.it;
            return *this;
        }
        
        inline ExpressionLink &operator*() { return **it; }
        inline ExpressionLink *operator->() { return *it; }
        
        inline iterator &operator++()
        {
            if (it != coll->end()) {
                ++it;
            }
            return *this;
        }
        
        inline iterator operator++(int)
        {
            iterator cop(coll, it);
            
            if (it != coll->end()) {
                ++it;
            }
            
            return cop;
        }

        inline iterator &operator--()
        {
            if (it != coll->begin()) {
                --it;
            }
            return *this;
        }
        
        inline iterator operator--(int)
        {
            iterator cop(coll, it);
            
            if (it != coll->begin()) {
                --it;
            }
            
            return cop;
        }

        inline iterator &operator+=(int n)
        {
            it += n;
            return *this;
        }
        
        inline iterator &operator-=(int n)
        {
            it -= n;
           return *this;
        }
        
        inline __w64 int operator-(const iterator &rhs) const
          { return it - rhs.it; }
        
        inline iterator operator+(int rhs) const
          { return iterator(coll, it + rhs); }
        
        inline iterator operator-(int rhs) const
          { return iterator(coll, it - rhs); }
        
        inline bool operator==(const iterator &rhs) const
          { return it == rhs.it; }

        inline bool operator!=(const iterator &rhs) const
          { return it != rhs.it; }

        inline bool operator<=(const iterator &rhs) const
          { return it <= rhs.it; }

        inline bool operator>=(const iterator &rhs) const
          { return it >= rhs.it; }

        inline bool operator<(const iterator &rhs) const
          { return it < rhs.it; }

        inline bool operator>(const iterator &rhs) const
          { return it > rhs.it; }
        
    private:
        std::vector<ExpressionLink *> *coll;
        std::vector<ExpressionLink *>::iterator it;
        
    private:
        friend class LinkCollection;
    };


    class name_iterator : public std::iterator_traits<ExpressionLink *>
    {
    public:
        name_iterator() : coll(0) {}
        
        name_iterator(std::vector<ExpressionLink *> *coll_, std::vector<ExpressionLink *>::iterator start, const std::string &name_)
            : coll(coll_), it(start), name(name_)
        {
            while (it != coll->end() && (*it)->name != name) {
                ++it;
            }
        }

        name_iterator &operator=(const name_iterator &i)
        {
            if (&i == this) {
                return *this;
            }
            
            coll = i.coll;
            it = i.it;
            name = i.name;
            return *this;
        }
        
        inline ExpressionLink &operator*() { return **it; }
        inline ExpressionLink *operator->() { return *it; }
        
        inline name_iterator &operator++()
        {
            if (it != coll->end()) {
                do {
                    ++it;
                } while (it != coll->end() && (*it)->name != name);
            }
            return *this;
        }
        
        inline name_iterator operator++(int)
        {
            name_iterator cop(coll, it, name);
            
            if (it != coll->end()) {
                do {
                    ++it;
                } while (it != coll->end() && (*it)->name != name);
            }
            
            return cop;
        }

        inline name_iterator &operator--()
        {
            if (it != coll->begin()) {
                do {
                    --it;
                } while (it != coll->begin() && (*it)->name != name);
            }
            return *this;
        }
        
        inline name_iterator operator--(int)
        {
            name_iterator cop(coll, it, name);
            
            if (it != coll->begin()) {
                do {
                    --it;
                } while (it != coll->begin() && (*it)->name != name);
            }
            
            return cop;
        }

        inline bool operator==(const name_iterator &rhs) const
          { return it == rhs.it; }

        inline bool operator!=(const name_iterator &rhs) const
          { return it != rhs.it; }

        inline bool operator<=(const name_iterator &rhs) const
          { return it <= rhs.it; }

        inline bool operator>=(const name_iterator &rhs) const
          { return it >= rhs.it; }

        inline bool operator<(const name_iterator &rhs) const
          { return it < rhs.it; }

        inline bool operator>(const name_iterator &rhs) const
          { return it > rhs.it; }
        
    private:
        std::vector<ExpressionLink *> *coll;
        std::vector<ExpressionLink *>::iterator it;
        std::string name;
        
    private:
        friend class LinkCollection;
    };

    
public:
    LinkCollection(std::vector<ExpressionLink *> &links_) : links(&links_) {}
    
    iterator begin() { return iterator(links, links->begin()); }
    iterator end()   { return iterator(links, links->end()); }

    name_iterator begin(const std::string &name) { return name_iterator(links, links->begin(), name); }
    name_iterator end(const std::string &name)   { return name_iterator(links, links->end(), name); }
    
    size_t size() const { return links->size(); }
    
    bool empty() const { return links->empty(); }

    iterator erase(iterator i) { return iterator(links, links->erase(i.it)); }
    
    ExpressionLink *operator[](size_t i) { return (*links)[i]; }
    
private:
    std::vector<ExpressionLink *> *links;
};


class ExpressionNode
{
public:
    ExpressionNode(RawStrokeGroup *strokes_) : strokes(strokes_), scope_depth(0) {}
    ~ExpressionNode()
    {
        while (!links.empty()) {
            delete links.back();
            links.pop_back();
        }
    }
    
    void AddTag(const std::string &tag);
    bool HasTag(const std::string &tag);
    void RemoveTag(const std::string &tag);
    
    void LinkTo(ExpressionNode *to, const std::string &name, double parm = 0.0);
    bool HasLink(ExpressionNode *to, const std::string &name);
    void RemoveLink(ExpressionNode *to, const std::string &name);
    LinkCollection GetLinks();

    const RawStrokeGroup &Strokes() const { return *strokes; }
    
    void AddMatch(const Match &match);
    std::vector<Match> &GetMatches() { return reco_results; }
    
    std::vector<std::list<unsigned> > &GetInputOrder() { return input_order; }
    void SetInputOrder(const std::vector<std::list<unsigned> > &io) { input_order = io; }
    void PushInputOrder(const std::list<unsigned> &io) { input_order.push_back(io); }
    
    unsigned scope_depth;
    
private:
    RawStrokeGroup *strokes;
    std::list<std::string> tags;
    std::vector<ExpressionLink *> links;
    
    std::vector<Match> reco_results;
    
    std::vector<std::list<unsigned> > input_order;
};

}


#endif

