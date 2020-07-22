enum flatten_action_overwrite {overwrite_image,return_correction_array};
enum flatten_action_xy {x_and_y,x_only,y_only};
enum flatten_action_linear {linear,nonlinear};
float *flatten_image(float *,int,int,
		     enum flatten_action_overwrite,
		     enum flatten_action_xy,
		     enum flatten_action_linear);
