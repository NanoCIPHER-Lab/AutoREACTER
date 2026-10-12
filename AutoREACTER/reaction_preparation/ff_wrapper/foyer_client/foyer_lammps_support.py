from pymatgen.io.lammps.data import LammpsData

class LammpsCombiner:

    @staticmethod
    def combine_lammps_files(file_list, output_filename= "force_field.data"):
        # Load the first data file as the base system
        combined_data = LammpsData.from_file(file_list[0])
    
        # Loop through and append the remaining files
        for next_file in file_list[1:]:
            additional_data = LammpsData.from_file(next_file)
            
            # combine() automatically merges atom types, shifts atom/bond IDs, 
            # and expands the simulation box to fit everything.
            combined_data = combined_data.combine(additional_data)
            
        # Write out the final combined LAMMPS data file
        combined_data.write_file(output_filename)
        print(f"Successfully combined {len(file_list)} files into: {output_filename}")

# # --- Example Usage ---
# # List the paths to the data files you want to combine
# files_to_merge = ["molecule1.data", "molecule2.data", "molecule3.data"]

# LammpsCombiner.combine_lammps_files(files_to_merge, "combined_system.data")