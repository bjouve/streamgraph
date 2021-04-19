The data files contain the whole data sets (6 days for the high school data and 2 days for the primary school) 

The results files contain 2 types of folders : the first one (input_W) contains the weighted graphs that constitute the input of the ItRich algorithm at each window. Each file is written like W_day_time (starting time of the window). Inside each file we can find an edge list (u,v,w) where u and v are the labels of the nodes and w their topological weight averaged during that window.  


The second folder (RC) contains the corresponding segmentation obtained by the ItRich Algorithm. For each time step, we have a list of rich clubs, ordered from the densest to the sparsest. The line indicating a rich club starts with "RC" and is followed by the labels of the nodes inside that rich club. This line is followed by the informations over the quality measure ( the line that starts with "Q") The corresponding numerical values are the average value and the standard deviation computed over a certain number of null models (N_null in the main.py) 
From these files, we can calculate the dense and the sparse part by setting a threshold over which the corresponding rich clubs are considered insides the dense part, and the rest of the nodes as the sparse part. 
The default value for the results showed in the article is 0.

Notes :

The user has to change the root folder (stored in the variable dir in the beginning of the script) in order to successfully execute the code. 

This code corresponds to a primary version and is far from being optimal in terms of time complexity.  