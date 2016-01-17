#pragma once

#include <thread>
using namespace std;

class ThreadGuard
{
	thread & t;
public:
	explicit ThreadGuard(thread& t_) : t(t_) {
	
	}
	~ThreadGuard();
	ThreadGuard(ThreadGuard const &) = delete;
	ThreadGuard& operator=(ThreadGuard const&) = delete;
};

