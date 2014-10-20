#include "log.h"


#define LOGFILE "C:\\libmathrec.log"


namespace scg
{


const log_stream::LogLevel log_stream::DefaultLevel(0);
const log_stream::LogLevel log_stream::InfoLevel(-5);
const log_stream::LogLevel log_stream::WarnLevel(3);
const log_stream::LogLevel log_stream::ErrorLevel(5);
const log_stream::LogLevel log_stream::CriticalLevel(10);
const std::string log_stream::endl("\n");

log_stream log;


}

