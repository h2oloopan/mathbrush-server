#ifndef LOG_H_
#define LOG_H_

#include <algorithm>
#include <fstream>


namespace scg
{


class log_stream : public std::ofstream
{
public:
    struct LogLevel {
        int level;
        
        LogLevel() : level(0) {}
        explicit LogLevel(int l) : level(l)
        {
            level = std::min(std::max(level, -10), 10);
        }
        
        bool operator<=(const LogLevel &rhs) const
        {
            return level <= rhs.level;
        }
    };
    
    static const LogLevel DefaultLevel;
    static const LogLevel InfoLevel;
    static const LogLevel WarnLevel;
    static const LogLevel ErrorLevel;
    static const LogLevel CriticalLevel;

    static const std::string endl;    
    
public:
    log_stream() : std::ofstream("C:\\libmathrec.log") {}

    template <typename T>
    log_stream &operator <<(T t)
    {
        if (thres <= level) {
            static_cast<std::ofstream &>(*this) << t;
        }
        return *this;
    }

    template <typename T>
    log_stream &operator <<(const T &t)
    {
        if (thres <= level) {
            static_cast<std::ofstream &>(*this) << t;
        }
        return *this;
    }
    
    log_stream &operator <<(const LogLevel &lev)
    {
        level = lev;
        return *this;
    }
    
    void threshold(const LogLevel &th)
    {
        thres = th;
    }
    
private:
    LogLevel level;
    LogLevel thres;
};


extern log_stream log;


}


#endif

