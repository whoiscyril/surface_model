#include "bulk.h"
#include "Input_Parser.h"
int main(int argc, char const *argv[])
{
	std::string filename = "input.in";
	//Creating the bulk model
	bulk_energy(get_input_coordinates(filename), get_input_species(filename), get_input_buckingham(filename));
	return 0;
}
