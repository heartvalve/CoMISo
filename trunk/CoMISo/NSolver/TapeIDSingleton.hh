/*
 * TapeIDSingleton.hpp
 *
 *  Created on: Jan 4, 2013
 *      Author: kremer
 */

#ifndef TAPEIDSINGLETON_HPP_
#define TAPEIDSINGLETON_HPP_

#include <string>

class TapeIDSingleton {
public:
    static TapeIDSingleton* Instance() {
        if(reference_ == NULL) {
            reference_ = new TapeIDSingleton();
        }
        return reference_;
    }

    short int uniqueTapeID() {
        short int id = tape_count_;
        ++tape_count_;
        return id;
    }
private:
    TapeIDSingleton() : tape_count_(0) {}
    TapeIDSingleton(const TapeIDSingleton&) {}
    ~TapeIDSingleton() {}
    static TapeIDSingleton* reference_;
    short int tape_count_;
};

#endif /* TAPEIDSINGLETON_HPP_ */
