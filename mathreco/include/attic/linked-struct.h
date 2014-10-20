#ifndef LINKED_STRUCT_H_
#define LINKED_STRUCT_H_


namespace scg
{


template <class T>
struct linked_struct
{
    linked_struct() : next(0) { }
    
    virtual ~linked_struct()
    {
        delete next;
    }
    
    T *next;
};


template <class T>
void
insert_linked_struct(T **head, T *entry)
{
    if (!*head) {
        *head = entry;
    }
    else {
        entry->next = *head;
        *head = entry;
    }
}


}


#endif
