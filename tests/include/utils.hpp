#include <sstream>

struct pvst_ss
{
    std::ostringstream  ss;
    std::string str() {return ss.str();}
    template<typename T>
    pvst_ss& operator<<(const T &t)
    {
        this->ss << t;
        return *this;
    }
};

#define STR(X) ( ( pvst_ss() << X ).str() )


uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return (((uint64_t)hi << 32) | lo) / 1000;
}
