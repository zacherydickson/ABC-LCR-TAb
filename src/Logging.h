#ifndef LOGGING_H
#define LOGGING_H

#pragma GCC diagnostic ignored "-Wformat-security"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>

namespace logger {

    enum LogLevel {ERROR = 1, WARNING = 2, INFO = 3, DEBUG = 4};

    extern int Verbosity;

    template<typename ... Args>
    void Log(const char * format, int lvl, Args ... args){
        if(lvl > Verbosity || lvl <= 0){
            return;
        }
        std::string header("[");
        //char buffer[100];
        //header += buffer;
        switch(lvl){
            case ERROR:
                header += "ERROR";
                break;
            case WARNING:
                header += "WARNING";
                break;
            case INFO:
                header += "INFO";
                break;
            case DEBUG:
            default:
                header += "DEBUG";
                break;
        }
        header += "] ";
        int size_s = std::snprintf( nullptr, 0, (header + format).c_str(), args ... ) + 1; // Extra space for '\0'
        if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
        auto size = static_cast<size_t>( size_s );
        std::unique_ptr<char[]> buf( new char[ size ] );
        std::snprintf( buf.get(), size, (header + format).c_str(), args ... );
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::stringstream stream;
        stream << std::put_time(&tm,"(%m-%d %H:%M:%S) ") << std::string(buf.get(), buf.get() + size - 1) << "\n";
        std::cerr << stream.str();
    }

}

#endif //LOGGING_H


