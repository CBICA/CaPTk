/** 
\file TrainingModuleResult.h

\brief This header file describes the trainingmodule result struct.
Basically just a boolean for success or failure, an error or success message, and an additional optional context field.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2021 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/

#pragma once
#include <string>

struct TrainingModuleResult
{
	bool success; // True if the operation succeeded, false otherwise
	std::string message; // Should be only what is directly relevant to what was just done
	std::string context; // Include caught exceptions, etc, here
};
