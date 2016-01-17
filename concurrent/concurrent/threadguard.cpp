#include "threadguard.h"


ThreadGuard::~ThreadGuard()
{
	if (t.joinable()) {
		t.join();
	}
}
