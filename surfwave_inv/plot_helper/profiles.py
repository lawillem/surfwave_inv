import numpy as np

def gen_depth_profiles(z_arr, property_arr, thk_arr):
    #z_arr, the depths at which you want output values
    #property_arr with size N contains the property we want to be interpolated
    #thk_arr is the layer thickness array. Either size N-1 or N with last value being np.inf, explicitly indicating infinite halfspace
    
    N = property_arr.size
    
    if thk_arr.size != N-1:
        if thk_arr.size != N or thk_arr[-1] != np.inf:
            raise Exception("Inconsistent thickness array provided")
        
        thk_arr = thk_arr[0:-1]
        
    #at this point thk_arr is size N-1 for sure.
    
    depth_interfaces = np.cumsum(thk_arr)
    
    ret = np.zeros_like(z_arr)
    for ilayer in xrange(N-1):
        
        bot_layer = depth_interfaces[ilayer]
        top_layer = bot_layer - thk_arr[ilayer]
        
        #find z within this window
        z_ind = np.logical_and(z_arr>=top_layer, z_arr<=bot_layer )
        ret[z_ind] = property_arr[ilayer]
        
    #take care of any points in z_arr in the bottom-most layer
    ret[z_arr > depth_interfaces[-1]] = property_arr[-1]
    return ret