import numpy as np 

myarr = []
degree_dist_arr = []
for i in range(3):
    myarr.append(np.arange(5))
    degree_dist_arr.append(np.arange(5))

def get_uniques(arr):
    b = []
    for a in arr:
        b.append(list(a))
    print(b)
    print(np.unique(b))

def mean_degree_dist(unique_degree_arr, degree_dist_arr):
    b = []
    for a in unique_degree_arr:
        b.append(list(a))
    uds = np.unique(b)
    ud_dict = {}
    deg_dict = {}
    av_degree_dist = []
    for ud in uds:
        ud_dict[ud] = 0
        deg_dict[ud] = 0

    for i,(ud_arr,deg_arr) in enumerate(zip(unique_degree_arr,degree_dist_arr)):
        for ud,deg_dist in zip(ud_arr, deg_arr):
            ud_dict[ud] = ud_dict[ud] + 1
            deg_dict[ud] = deg_dict[ud] + deg_dist

    for ud in ud_dict.keys() :
        print(f"Average degree dist[{ud}] = {deg_dict[ud]}/{ud_dict[ud]}")
        av_degree_dist.append(deg_dict[ud]/ud_dict[ud])
    print(f"Unique degree array: {unique_degree_arr}")
    print(f"Degree dist array: {degree_dist_arr}")
    print(f"ud_dict: {ud_dict}\ndeg_dict:{deg_dict}")
    print(f"deg dict values: {deg_dict.values()}")

    print(f"Average degrees: {av_degree_dist}")
    return av_degree_dist

    

    # print()





# print(myarr)
# # print(myarr[:][0])
# print(f"Outside:{np.unique(myarr)}")
# get_uniques(myarr)
# print(myarr[-1] + [0,0,0,0,-1])
# degree_dist_arr[-1] = degree_dist_arr[-1] + [0,0,0,0,-2]
degree_dist_arr[-1] = [0,0]
mean_degree_dist(myarr, degree_dist_arr)





