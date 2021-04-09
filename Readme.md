# TilingSquare
Program to discard all the possible 4,5-sides tilings of the square with seven (maybe more) pieces.

## INTRODUCTION
\
Let **S** be a square and **T** be a convex polygon. We say that **S** can be tiled by **n** copies of **T** if there are **n** convex polygons, all congruent to **T** such that the union of the copies are equal to **S** and the the polygons have disjoint interiors.

If **n** is an odd number, then it is easy to see that **S** can be split into **n** congruent rectangles. It is not known in general if **S** can be tiled by **n** copies of **T** when **T** is not a rectangle. The standing conjecture is as follows.

**Conjecture**. If **n >= 3** is an odd integer, then a square **S** may be tiled by **n** congruent copies of a convex polygon **T** only if **T** is a rectangle.

To see references and more details about this problem and the idea which lead to code this see: **aqui va la liga al articulo**

## What does this program?

### Pre-preparation of data
\
First we need the files with the graphs associated with the possible admissible tilings. We previously provided this files, the ones without extension and name **XX**pm**Y**c**Z**, and use it as input depending of in how many tiles we want to split the square.

(**XX** = number of tiles plus five)

Example for 7 pieces:
```python
Graphs7 = load_plantri(filename='12pm4c3')
```

## Usage
### One notebook for each amount of tiles
\

In each file, named TilingSearch**N**.ipynb, we first add a new graph for each possible distinguished vertex.

```python
inputs = tqdm(Graphs7)

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(add_distinguished)(G) for G in inputs)

Graphs7 = [G for DGraph in processed_list for G in DGraph]

print(f'Constructed {len(Graphs7)} graphs with a distinguished vertex.')
save_data(Graphs7,'all_tiling_grahps_7.txt')
```

After that the remaining graphs pass through some filters based in some configurations that is known are not possible.

```python
#Graphs7 = load_data('all_tiling_grahps_7.txt')

inputs = tqdm(Graphs7)

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(first_filter)(G) for G in inputs)

Graphs7 = [Graphs7[i] for i in range(len(Graphs7)) if processed_list[i]]

print(f'Left with {len(Graphs7)} graphs after first filter.')
save_data(Graphs7,'filtered_tiling_grahps_7.txt')
```

Notice that we have commented the line in which it can be loaded the previous results instead of calculate everything from the start.
\
\
The last part consists in exhaust the remaining graphs using a deep-search algorithm to assign angles and sides to each piece and see if it could be an achievable tiling.

```python
#Graphs7_4 = load_data('quadrilateral_tiling_grahps_7.txt')

save_ang_perms(4,'4_perms.txt')

inputs = tqdm(Graphs7_4)
numSides = 4

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(search)(ig,G,numSides) for ig,G in enumerate(inputs))

FinalNodes = [Nodes for Nodes in processed_list if len(Nodes)>0]

print(f'Left with {len(FinalNodes)} graphs.')
```

## Usage

You only need to run the notebook in order. If you have run a previous filter, you can uncomment the first line of the next box to use that calculated data.

## Output
The output is a list of lists. Each one containing objects that says which graph it is (with the indices in the input), with which permutation of angle-types (with the index it has in the generated perms file), and some additional information used in the process.
\
\
In the case the output is empty, means that all the cases were exhausted.
