#include <stdlib.h>
#include <unistd.h>
#include <iostream>

using namespace std;

int main (int argc, char *argv[]) { 
	setenv("LD_LIBRARY_PATH", "../lib", 1);
	pid_t pid = fork();  //  creates a second process, an exact copy of the current one
	if (pid==0)  {  // this is exectued in the child process
	    if (execvp("./CaPTk", argv))   // execvp() returns only if lauch failed
		 cout << "Couldn't run CaPTk" << endl;   
	}
	else {  // this is executed in the parent process 
	    if (pid==-1)   //oops ! This can hapen as well :-/
		 cout << "Process launch failed";  
	    else cout << "I launched process "<<pid<<endl; 
	}
} 
