# WHAT IS IT?

This is an implementation of the algorithms proposed in the article [Dissecting the square into seven or nine congruent parts](https://www.sciencedirect.com/science/article/pii/S0012365X22000061?casa_token=oJryTNyL9YcAAAAA:sPaHFUhVSasikQ4vtMCgTpSS4rrfBJDjAdFnFD-nWQ0KUFSbmndg0zVjQ-P8UcIYNlEJpn6xr3ZX). We recommend reading it first.

The executables are provided on **Jupyter** notebooks, these are divided into **search** and **proof** files. In the first type of file, the main process is executed. In the second, the user can see a step by step explanation of what the deep search does for a given graph. To provide a more detailed explanation of what these files do, we use the square case as an example. The notebooks for other variants (tiling a rectangle or using equiangular pieces) are similar but use filters that fit their conditions. It is important to note that **joblib.Parallel** is used to parallelize the process.

In the output, the current graph is drawn (not necessarily a planar drawing) and the proof follows.

## SEARCH FILES

The search files begin by loading the graphs generated by **plantri**:
``` python
Graphs7 = load_plantri(filename='plantri/12pm4c3')
print(f'Loaded {len(Graphs7)} graphs from plantri.')
```
The **plantri** files are not included in this repository since some of these are quite large. To see more about plantri and how to use it see the [documentation](https://users.cecs.anu.edu.au/~bdm/plantri/plantri-guide.txt).
\
The notebook continues creating a copy of each graph for each possible selection of the distinguished vertices (see remark in section 3 of the paper): 
``` python
inputs = tqdm(Graphs7)
# add_distinguished is executed in each graph of Graphs7
if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(add_distinguished)(G) for G in inputs)
Graphs7 = [G for DGraph in processed_list for G in DGraph]
print(f'Constructed {len(Graphs7)} graphs with a distinguished vertex.')
```
After that, we check that condition 2 of Lemma 5 is fulfilled:
``` python
inputs = tqdm(Graphs7)
if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(first_filter)(G) for G in inputs)
Graphs7 = [Graphs7[i] for i in range(len(Graphs7)) if processed_list[i]]
print(f'There are {len(Graphs7)} graphs left after first filter.')
```
For the remaining graphs, the cases are separated by the possible number of sides of **T**:
``` python
inputs = tqdm(Graphs7)
if __name__ == "__main__":
    #here we are separating the graphs that could represent a tiling with a polygon with five sides.
    processed_list = Parallel(n_jobs=num_cores)(delayed(has_min_deg)(G,5) for G in inputs)
Graphs7_5 = [Graphs7[i] for i in range(len(Graphs7)) if processed_list[i]]
print(f'There are {len(Graphs7_5)} pentagonal graphs.')
```
At this point, we have a list of remaining possible graphs for each number of sides. 
For quadrilateral tiles, we apply a second filter:

``` python
inputs = tqdm(Graphs7)
if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(second_filter)(G) for G in inputs)
Graphs7_4 = [Graphs7[i] for i in range(len(Graphs7)) if processed_list[i]]
print(f'There are {len(Graphs7_4)} quadrilateral graphs left after second filter.')
save_data(Graphs7_4,'graphs/filtered_7_4.txt')
```
This is the function that checks that Lemma 5 condition 3 is fulfilled.
\
At the end, the deep search algorithm is applied over each remaining graph:
``` python
Graphs7_4 = load_data('filtered_7_4.txt')
#Here the possible labelings are generated
save_ang_perms(4,'4_perms.txt')
inputs = tqdm(Graphs7_4)
if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(search)(ig,G,'4_perms.txt') for ig,G in enumerate(inputs))
FinalNodes = [Nodes for Nodes in processed_list if len(Nodes)>0]
print(f'There are {len(FinalNodes)} graphs left after performing a deep search for each quadrilateral graph.')
```

## PROOF FILES

As mentioned before, you can generate a step by step explanation of what the deep search does for a given graph. These files start once the preliminary filters have been applied just as in the search files. After that, a proof for a given graph is generated:
``` python
import random
save_ang_perms(4,'4_perms.txt')
Graphs5_4 = load_data('graphs/filtered_5_4.txt')
#Choose a random graph
ig=random.randrange(len(Graphs5_4))
print('Incidence graph',ig)
G=Graphs5_4[ig] 
G.draw()
search(ig,G,'4_perms.txt',PrintProof=True)
```
In the output, the current graph is drawn (not necessarily a planar drawing) and the proof follows. We recommend seeing a previously generated proof in one of these files. When a graph has possible equations to be realized, these equations are shown.

# NOTES ABOUT RESULTS

## EQUIANGULAR CASE

In the paper, we mention that some graphs are not filtered by the algorithm but are unrealizable. Every graph that is not filtered can be found in the search files. We list the graphs that do not produce a valid tiling here, using its respective number assigned in the search file:

<ul>
    <li>5 tiles, 3 sides: 510, 1298, 1586</li>
    <li>5 tiles, 4 sides: 6, 11, 60, 68, 72, 89, 91</li>
</ul> 
