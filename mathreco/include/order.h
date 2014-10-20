#ifndef ORDER_H_
#define ORDER_H_


#include "reco-types.h"

#define NUM_ORDERS 3


namespace scg
{


extern const unsigned X_ORDER;
extern const unsigned Y_ORDER;
extern const unsigned TIME_ORDER;

typedef bool (*order_fn)(const segment *, const segment *);

extern order_fn orders[NUM_ORDERS];


}


#endif
