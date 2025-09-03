#include "common.h"
#include "LaserGenerator.h"
#include "MuGenerator.h"
#include "RootManager.h"
#include "RunManager.h"

// Meyer's Singleton Pattern
RunManager *RM_ptr = &RunManager::GetInstance();
LaserGenerator *lsr_ptr = &LaserGenerator::GetInstance();
MuGenerator *Mu_ptr = &MuGenerator::GetInstance();
RootManager *ROOT_ptr = &RootManager::GetInstance();


int main(int argc, char* argv[]) {
    // Check if a macro file is provided
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <macro_file_path>" << std::endl;
        return 1;
    }

    // Get the macro file path from the command line arguments
    std::string macro_file_path = argv[1];

    RM_ptr->ReadCommandFile(macro_file_path);

    RM_ptr->SolveOBE();

    std::cout << "\nSimulation completed." << std::endl;
    return 0;
}
