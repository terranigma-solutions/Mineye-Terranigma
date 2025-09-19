def plutonic_perturbations(level_of_sampling_points, level_of_sampling_orientations, 
                       section_dict,
                       plutonite_orientations_df, plutonite_points_df,
                       topography_path,
                       results_folder):

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
    
    methods_of_assigning_dip_angles = {"purely random": 1, "mildly random": 2, "little random": 3} # controlling randomness of dip angles

    for sampling in level_of_sampling_points:
        for orientations in level_of_sampling_orientations:
            for method_name, method_value in methods_of_assigning_dip_angles.items():

                # Create a descriptive string for the current parameter combination
                param_combination = f"Sampling Points: {sampling}, Sampling Orientations: {orientations}, Method: {method_name} (value: {method_value})"

                # # Check in the target folder if the image already exists
                # target_folder = os.path.join(main_results_folder, 
                #                         f"Sampling_Points_{int(sampling*100)}pct",
                #                         f"Sampling_Orientations_{int(orientations*100)}pct")
                # image_filename = os.path.join(target_folder, 
                #                             f"plutonite_model_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name.replace(' ', '_')}.png")
                # if os.path.exists(image_filename):
                #     print(f"Image already exists for {param_combination}, skipping computation.")
                #     continue  # Skip to the next iteration if image already exists


                # # Make a copy of the orientations DataFrame
                # plutonite_orientations_df = plutonite_orientations_df.copy()
                
                # """
                # As it is described now in the distribution of dip angles:
                # we assume plutonites are thicker in the West than in the East.
                # """

                # # Assign dip angles based on the chosen method
                # if method_name == "purely random":
                #     # Very random
                #     std = 5
                #     plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(80, std, size=27)     # indices 0 to 26
                #     plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(60, std, size=14)     # indices 27 to 40
                #     plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(60, std, size=2)      # indices 41 to 42
                #     plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(40, std, size=7)      # indices 43 to 49
                #     plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(60, std, size=1)      # index 50 only
                # elif method_name == "mildly random":
                #     # Midly random
                #     std = 2.5
                #     plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(80, std, size=27)     # indices 0 to 26
                #     plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(60, std, size=14)     # indices 27 to 40
                #     plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(60, std, size=2)      # indices 41 to 42
                #     plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(40, std, size=7)      # indices 43 to 49
                #     plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(60, std, size=1)      # index 50 only
                # elif method_name == "little random":
                #     # Little random
                #     std = 1
                #     plutonite_orientations_df.loc[0:26,  'dip']    = np.random.normal(80, std, size=27)     # indices 0 to 26
                #     plutonite_orientations_df.loc[27:40, 'dip']    = np.random.normal(60, std, size=14)     # indices 27 to 40
                #     plutonite_orientations_df.loc[41:42, 'dip']    = np.random.normal(60, std, size=2)      # indices 41 to 42
                #     plutonite_orientations_df.loc[43:49, 'dip']    = np.random.normal(40, std, size=7)      # indices 43 to 49
                #     plutonite_orientations_df.loc[50,    'dip']    = np.random.normal(60, std, size=1)      # index 50 only

                
                # # Clip the values to ensure they're between 0 and 90
                # plutonite_orientations_df['dip'] = np.clip(plutonite_orientations_df['dip'], 0, 90)


    #             # Calculate number of data points to sample
    #             num_points_to_sample_points = int(len(plutonite_orientations_df) * sampling)
    #             num_points_to_sample_orientations = int(len(plutonite_points_df) * orientations)

    #             # Sample the data based on the current sampling level without random state set (!)
    #             sampled_orientations_df = plutonite_orientations_df.sample(n=num_points_to_sample_points).reset_index(drop=True)
    #             sampled_points_df = plutonite_points_df.sample(n=num_points_to_sample_orientations).reset_index(drop=True)
    #             print(f"Processing: Sampling level (points) = {sampling}, Sampling level (orientations) = {orientations}, Method = {method_name} (value: {method_value})")

    #             # Save the sampled data to temporary CSV files
    #             sampled_orientations_path = "temp_sampled_orientationinputpoints.csv"
    #             sampled_points_path = "temp_sampled_formationinputpoints.csv"
    #             sampled_orientations_df.to_csv(sampled_orientations_path, index=False)
    #             sampled_points_df.to_csv(sampled_points_path, index=False)

    #             # Create geological model
    #             plutonite_geo_model = gp.create_geomodel(
    #                 project_name='AOI',
    #                 extent=extent,
    #                 # refinement=6,
    #                 resolution=resolution,
    #                 importer_helper=gp.data.ImporterHelper(
    #                     path_to_orientations=os.path.join(BASE_DIR,   "temp_sampled_orientationinputpoints.csv"),
    #                     path_to_surface_points=os.path.join(BASE_DIR, "temp_sampled_formationinputpoints.csv"),
    #                 )
    #             )
                
    #             # Set layers
    #             gp.map_stack_to_surfaces(
    #                 gempy_model=plutonite_geo_model,
    #                 mapping_object={
    #                     "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    #                 }
    #             )

    #             # Set topography
    #             gp.set_topography_from_file(grid=plutonite_geo_model.grid, filepath=path_to_topography_reduced)

    #             # Set cross sections
    #             gp.set_section_grid(plutonite_geo_model.grid, section_dict)
                
    #             # Compute plutonite model
    #             plutonite_model = gp.compute_model(plutonite_geo_model)

    #             # Define the target subfolder for this iteration
    #             target_folder = os.path.join(main_results_folder, 
    #                                     f"Sampling_Points_{int(sampling*100)}pct",
    #                                     f"Sampling_Orientations_{int(orientations*100)}pct")

    #             # Create image filename with path to target folder
    #             image_filename = os.path.join(target_folder, 
    #                                         f"plutonite_model_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name.replace(' ', '_')}.png")
                
    #             # Create image filenames for cross sections
    #             W_E_cross_section_filename = os.path.join(target_folder, 
    #                                         f"W_E_plutonite_model_cross_sections_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name.replace(' ', '_')}.png")
    #             N_S_cross_section_filename = os.path.join(target_folder, 
    #                                         f"N_S_plutonite_model_cross_sections_sampling_points_{int(sampling*100)}pct_sampling_orientations_{int(orientations*100)}pct_method_{method_name.replace(' ', '_')}.png")

    #             # Visualise results from above and display parameter combination
    #             plutonite_lith_block = plutonite_geo_model.solutions.raw_arrays.lith_block
    #             plutonite_lith_block_reshaped = plutonite_lith_block.reshape(64,64,64)
    #             plutonite_id = 1
    #             plutonite_mask = plutonite_lith_block_reshaped == plutonite_id
    #             plutonite_count_2d = np.sum(plutonite_mask, axis=2)
    #             plutonite_count_2d_masked = np.ma.masked_where(plutonite_count_2d == 0, plutonite_count_2d)
                
    #             # Create the 2D plot using imshow with masked zeros
    #             fig, ax = plt.subplots(figsize=(12, 10))

    #             # Use imshow for a pixel-based representation with masked zeros
    #             im = ax.imshow(plutonite_count_2d_masked.T, 
    #                         extent=[min_x/1000, max_x/1000, min_y/1000, max_y/1000], 
    #                         cmap='plasma', 
    #                         origin='lower', 
    #                         interpolation='nearest', 
    #                         alpha=0.9)

    #             # Add colorbar
    #             cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    #             cbar.set_label('Number of Plutonite Voxels in Vertical Stack', fontsize=12)

    #             # Set labels and title
    #             ax.set_xlabel('X Coordinate (km)', fontsize=12)
    #             ax.set_ylabel('Y Coordinate (km)', fontsize=12)
    #             ax.set_title(f'Plutonite Thickness\n{param_combination}', fontsize=14, pad=20)

    #             # Set aspect ratio to be equal
    #             ax.set_aspect('equal')

    #             # Add grid
    #             ax.grid(True, alpha=0.3, linestyle='--', color='black')

    #             # draw cross section lines
    #             for section_name, (start, end, res) in section_dict.items():
    #                 ax.plot([start[0]/1000, end[0]/1000], [start[1]/1000, end[1]/1000], 
    #                         color='red', linestyle='--', linewidth=1.5, alpha=0.8)
    #                 mid_x = (start[0] + end[0]) / 2000
    #                 mid_y = (start[1] + end[1]) / 2000
    #                 ax.text(mid_x, mid_y, section_name, color='white', fontsize=7, ha='center', va='center', backgroundcolor='black', alpha=0.7)
                
    #             plt.tight_layout()

    #             # Save the image to the appropriate subfolder
    #             fig.savefig(image_filename, dpi=300, bbox_inches='tight')
    #             plt.close()
    #             print(f"Saved image: {image_filename}")

    #             # Create cross sections
    #             gpv.plot_2d(plutonite_geo_model, section_names=['W-E_1','W-E_2','W-E_3','W-E_4'], ve=15)
    #             # save the image
    #             plt.savefig(W_E_cross_section_filename, dpi=300, bbox_inches='tight')
    #             plt.close()

    #             gpv.plot_2d(plutonite_geo_model, section_names=['N-S_1','N-S_2','N-S_3','N-S_4'], ve=15)
    #             # save the image
    #             plt.savefig(N_S_cross_section_filename, dpi=300, bbox_inches='tight')
    #             plt.close()
                
    #             # Delete temporary CSV files
    #             if os.path.exists(sampled_orientations_path):
    #                 os.remove(sampled_orientations_path)
    #             if os.path.exists(sampled_points_path):
    #                 os.remove(sampled_points_path)

    #             # Display progress
    #             model_number = (list(level_of_sampling_points).index(sampling) * len(level_of_sampling_orientations) * len(methods_of_assigning_dip_angles)) + (list(level_of_sampling_orientations).index(orientations) * len(methods_of_assigning_dip_angles)) + list(methods_of_assigning_dip_angles.keys()).index(method_name) + 1
    #             total_models = len(level_of_sampling_points) * len(level_of_sampling_orientations) * len(methods_of_assigning_dip_angles)
    #             print(f"Completed model number {model_number} of {total_models}")

    # gpv.plot_section_traces(plutonite_geo_model)
    # print("\nAll perturbation results have been saved to organized subfolders!")
    print("Function plutonic_perturbations executed (function body commented out).")