#ifndef HELPER_FUNCS_H_INCLUDED
#define HELPER_FUNCS_H_INCLUDED 

#include <iostream> 

void input() { // Python style input function, useful for debug
    std::string buf;
    std::getline(std::cin, buf);
}
#endif 