import pandas as pd
import numpy as np
import os
import gempy as gp
import gempy_viewer as gpv
import matplotlib.pyplot as plt

def plutonic_perturbations(level_of_sampling_points, level_of_sampling_orientations, 
                       section_dict,
                       plutonite_orientations_df_path, plutonite_points_df_path,
                       topography_path,
                       main_results_folder,
                       dip_angle_mean, dip_angle_std, dip_angle_bounds,
                       extent, resolution):

    """
    This function creates multiple geological models of plutonites by varying the sampling levels of the
    contact points, orientations points, as well as providing three levels of randomness in assigning dip angles
    of the orientation points. These three levels are automatically defined in the function.
    The results are saved in an organized folder structure for easy access and comparison with images of the models
    in plan view and cross sections
    Parameters:
    - level_of_sampling_points: List or array of sampling levels for contact points (e.g., [0.9, 0.8, 0.7])
    - level_of_sampling_orientations: List or array of sampling levels for orientation points (e.g., [0.9, 0.8, 0.7])
    - section_dict: Dictionary defining cross sections with their start and end coordinates and resolution
    Returns:
    - None (saves the images and stores the geological models in organized folders)
    """

    # Create subfolder structure for all combinations to store the results later
    for sampling in level_of_sampling_points:
        sampling_folder = os.path.join(main_results_folder, f"Sampling_Points_{int(sampling*100)}pct")
        if not os.path.exists(sampling_folder):
            os.makedirs(sampling_folder)
        
        for orientations in level_of_sampling_orientations:
            orientations_folder = os.path.join(sampling_folder, f"Sampling_Orientations_{int(orientations*100)}pct")
            if not os.path.exists(orientations_folder):
                os.makedirs(orientations_folder)

    methods_of_assigning_dip_angles = {"Very_random": 1, "Mildly_random": 2, "Little_random": 3} # controlling randomness of dip angles

    for sampling in level_of_sampling_points:
        for orientations in level_of_sampling_orientations:
            for method_name, method_value in methods_of_assigning_dip_angles.items():

                # Create a descriptive string for the current parameter combination
                param_combination = f"Sampling Points: {sampling}, Sampling Orientations: {orientations}, Method: {method_name} (value: {method_value})"


                # Check in the target folder if the results already exists
                target_folder = os.path.join(main_results_folder, 
                                        f"Sampling_Points_{int(sampling*100)}pct",
                                        f"Sampling_Orientations_{int(orientations*100)}pct")
                image_filename = os.path.join(target_folder, 
                                            f"plutonite_model_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name}.png")
                if os.path.exists(image_filename):
                    print(f"Image already exists for {param_combination}, skipping computation.")
                    continue  # Skip to the next iteration if results already exists

                # Ensure the target folder exists
                os.makedirs(target_folder, exist_ok=True)

               # Create cross-section filenames with proper path
                W_E_cross_section_filename = os.path.join(target_folder, 
                    f"W_E_plutonite_model_cross_sections_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name}.png")

                N_S_cross_section_filename = os.path.join(target_folder, 
                    f"N_S_plutonite_model_cross_sections_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name}.png")

                # Ensure directory exists for cross-section files
                os.makedirs(os.path.dirname(W_E_cross_section_filename), exist_ok=True)
                os.makedirs(os.path.dirname(N_S_cross_section_filename), exist_ok=True)

                # Create image filename for this iteration
                image_filename = os.path.join(target_folder, 
                                            f"plutonite_model_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name}.png")

                # Make a copy of the orientations and points DataFrames
                plutonite_orientations_df = pd.read_csv(plutonite_orientations_df_path)
                plutonite_orientations_df = plutonite_orientations_df.copy()

                plutonite_points_df = pd.read_csv(plutonite_points_df_path)
                plutonite_points_df = plutonite_points_df.copy()


                # Assign dip angles based on the method
                if method_name == "Very_random":
                    # Very random
                    std = dip_angle_std["Very_random"]
                    plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(dip_angle_mean["West Plutonite"],       std, size=27)
                    plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=14)
                    plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=2)
                    plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(dip_angle_mean["South Plutonite"],      std, size=7)
                    plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=1)
                elif method_name == "Mildly_random":
                    # Mildly random
                    std = dip_angle_std["Mildly_random"]
                    plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(dip_angle_mean["West Plutonite"],       std, size=27)
                    plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=14)
                    plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=2)
                    plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(dip_angle_mean["South Plutonite"],      std, size=7)
                    plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=1)
                elif method_name == "Little_random":
                    # Little random
                    std = dip_angle_std["Little_random"]
                    plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(dip_angle_mean["West Plutonite"],       std, size=27)
                    plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=14)
                    plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=2)
                    plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(dip_angle_mean["South Plutonite"],      std, size=7)
                    plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(dip_angle_mean["North-East Plutonite"], std, size=1)

                
                # Clip the values to ensure they're between appropriate bounds
                plutonite_orientations_df['dip'] = np.clip(plutonite_orientations_df['dip'], dip_angle_bounds['Lower'], dip_angle_bounds['Upper'])


                # Calculate number of data points to sample
                num_points_to_sample_points = int(len(plutonite_orientations_df) * (1-sampling))
                num_points_to_sample_orientations = int(len(plutonite_points_df) * (1-orientations))

                # Randomly sample the data based on the current sampling level without random state set
                sampled_orientations_df = plutonite_orientations_df.sample(n=num_points_to_sample_points).reset_index(drop=True)
                sampled_points_df       = plutonite_points_df.sample(n=num_points_to_sample_orientations).reset_index(drop=True)
                
                # Display the current parameter combination
                print(f"Processing: Sampling level (points) = {sampling}, Sampling level (orientations) = {orientations}, Method = {method_name} (value: {method_value})")


                # Save the sampled data to temporary CSV files
                sampled_orientations_path = "temp_sampled_orientationinputpoints.csv"
                sampled_points_path       = "temp_sampled_formationinputpoints.csv"
                sampled_orientations_df.to_csv(sampled_orientations_path, index=False)
                sampled_points_df.to_csv(sampled_points_path, index=False)


                # Create geological model
                plutonite_geo_model = gp.create_geomodel(
                    project_name='Perturbed_Plutonite_Model',
                    extent=extent,
                    # refinement=6,
                    resolution=resolution,
                    importer_helper=gp.data.ImporterHelper(
                        path_to_orientations="temp_sampled_orientationinputpoints.csv",
                        path_to_surface_points="temp_sampled_formationinputpoints.csv",
                    )
                )
                
                # Set layers
                gp.map_stack_to_surfaces(
                    gempy_model=plutonite_geo_model,
                    mapping_object={
                        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
                    }
                )

                # Set topography
                gp.set_topography_from_file(grid=plutonite_geo_model.grid, filepath=topography_path)

                # Set cross sections
                gp.set_section_grid(plutonite_geo_model.grid, section_dict)
                
                # Compute plutonite model
                plutonite_model = gp.compute_model(plutonite_geo_model)


                # Extract the lithology voxels
                plutonite_lith_block          = plutonite_geo_model.solutions.raw_arrays.lith_block
                plutonite_lith_block_reshaped = plutonite_lith_block.reshape(64,64,64)

                # Locate which voxels contain plutonite
                plutonite_id = 1
                plutonite_mask = plutonite_lith_block_reshaped == plutonite_id

                # Calculate the number of plutonite voxels in each vertical stack (x,y)
                # and mask zeros for better visualization
                plutonite_count_2d = np.sum(plutonite_mask, axis=2)
                plutonite_count_2d_masked = np.ma.masked_where(plutonite_count_2d == 0, plutonite_count_2d)
                

                # Visualise the shape and thickness of plutonite blobs
                fig, ax = plt.subplots(figsize=(12, 10))

                min_x = extent[0]
                max_x = extent[1]
                min_y = extent[3]
                max_y = extent[2]

                im = ax.imshow(plutonite_count_2d_masked.T,
                            extent=[min_x/1000, max_x/1000, min_y/1000, max_y/1000], 
                            cmap='plasma', 
                            origin='lower', 
                            interpolation='nearest', 
                            alpha=0.9)

                cbar = plt.colorbar(im, ax=ax, shrink=0.8)
                cbar.set_label('Number of Plutonite Voxels in Vertical Stack', fontsize=12)

                ax.set_xlabel('X Coordinate (km)', fontsize=12)
                ax.set_ylabel('Y Coordinate (km)', fontsize=12)
                ax.set_title(f'Plutonite Thickness\n{param_combination}', fontsize=14, pad=20)

                ax.set_aspect('equal')

                ax.grid(True, alpha=0.3, linestyle='--', color='black')

                # Visualise where the crosssections are
                for section_name, (start, end, res) in section_dict.items():
                    ax.plot([start[0]/1000, end[0]/1000], [start[1]/1000, end[1]/1000], 
                            color='red', linestyle='--', linewidth=1.5, alpha=0.8)
                    mid_x = (start[0] + end[0]) / 2000
                    mid_y = (start[1] + end[1]) / 2000
                    ax.text(mid_x, mid_y, section_name, color='white', fontsize=7, ha='center', va='center', backgroundcolor='black', alpha=0.7)
                
                plt.tight_layout()

                # Save the image to the appropriate subfolder
                fig.savefig(image_filename, dpi=300, bbox_inches='tight')
                plt.close()


                # # Create and store cross sections in the appropriate subfolder
                """
                There is a problem with saving the crossections. I've been trying to fix it, but I
                haven't found the solution yet.
                """
                # # Add this line before saving cross-sections:
                # os.makedirs(target_folder, exist_ok=True)
                # # Create and save W-E cross sections
                # gpv.plot_2d(plutonite_geo_model, section_names=['W-E_1','W-E_2','W-E_3','W-E_4'], ve=15)
                # plt.savefig(W_E_cross_section_filename, dpi=300, bbox_inches='tight')
                # plt.close()

                # # Create and save N-S cross sections  
                # gpv.plot_2d(plutonite_geo_model, section_names=['N-S_1','N-S_2','N-S_3','N-S_4'], ve=15)
                # plt.savefig(N_S_cross_section_filename, dpi=300, bbox_inches='tight')
                # plt.close()

                

                # Delete temporary CSV files
                if os.path.exists(sampled_orientations_path):
                    os.remove(sampled_orientations_path)
                if os.path.exists(sampled_points_path):
                    os.remove(sampled_points_path)


# DO FORWARD MODELLING AND CALCULATE MISFIT GEOPHYSICS AND CALCULATE MISFIT GEOLOGY

# SAVE THE MODEL


                # Display progress
                model_number = (list(level_of_sampling_points).index(sampling) * len(level_of_sampling_orientations) * len(methods_of_assigning_dip_angles)) + (list(level_of_sampling_orientations).index(orientations) * len(methods_of_assigning_dip_angles)) + list(methods_of_assigning_dip_angles.keys()).index(method_name) + 1
                total_models = len(level_of_sampling_points) * len(level_of_sampling_orientations) * len(methods_of_assigning_dip_angles)
                
                print(f"Completed model number {model_number} of {total_models}")

    # Only plot section traces if a model was actually computed
    if 'plutonite_geo_model' in locals():
        gpv.plot_section_traces(plutonite_geo_model)
    
    print("\nAll perturbation results have been saved to organized subfolders!")