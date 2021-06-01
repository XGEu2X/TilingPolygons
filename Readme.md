# TilingPolygons
Program to list all the possible congruent pieces tilings of a Quadrilateral. Also includes a program to list all the possible equiangular pieces tilings of a square.

## INTRODUCTION
\
Let **S** be a square and **T** be a convex polygon. We say that **S** can be tiled by **n** copies of **T** if there are **n** convex polygons, all congruent to **T** such that the union of the copies are equal to **S** and the the polygons have disjoint interiors.

If **n** is an odd number, then it is easy to see that **S** can be split into **n** congruent rectangles. It is not known in general if **S** can be tiled by **n** copies of **T** when **T** is not a rectangle. The standing conjecture is as follows.

**Conjecture**. If **n >= 3** is an odd integer, then a square **S** may be tiled by **n** congruent copies of a convex polygon **T** only if **T** is a rectangle.

With this code, we could check positively the conjecture for **n** equal to **7** and **9**. Also the program works for other variants of the conjecture.

To see references and more details about this problem and the idea which lead to code this see: https://arxiv.org/abs/2104.04940

## What does this program?

### Pre-preparation of data
\
First we need the files with the graphs associated with the possible admissible tilings. The user can generate this graphs using **plantri** (http://users.cecs.anu.edu.au/~bdm/plantri/), making all the graphs with **n+5** vertices and puting conditions on minimum degree and conexity depending of the type of Tiling.

The file names of the must contain the next format: **XX**pm**Y**c**Z**. Where:

**XX** = number of tiles plus five
**Y** = minimum degree (if it is requiered)
**Z** = minimum conexity number

The code reads this files to work, and are absolutely necessarily.

Example for 7 congruent pieces:
```python
Graphs7 = load_plantri(filename='12pm4c3')
```

## Usage

### One notebook for each amount of tiles

In each file, named TilingSearch**N**.ipynb, we first add a new graph for each possible distinguished vertex.

```python
inputs = tqdm(Graphs7)

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(add_distinguished)(G) for G in inputs)

Graphs7 = [G for DGraph in processed_list for G in DGraph]

print(f'Constructed {len(Graphs7)} graphs with a distinguished vertex.')
```

After that the remaining graphs pass through some filters based in some configurations that is known are not possible.

```python
inputs = tqdm(Graphs7)

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(first_filter)(G) for G in inputs)

Graphs7 = [Graphs7[i] for i in range(len(Graphs7)) if processed_list[i]]

print(f'Left with {len(Graphs7)} graphs after first filter.')
```

The last part consists in exhaustively check the remaining graphs using a deep-search algorithm to assign angles and sides to each piece and see if it could be an bildable tiling.

```python
Graphs7_4 = load_data('filtered_7_4.txt')

save_ang_perms(4,'4_perms.txt')

inputs = tqdm(Graphs7_4)

if __name__ == "__main__":
    processed_list = Parallel(n_jobs=num_cores)(delayed(search)(ig,G,'4_perms.txt') for ig,G in enumerate(inputs))

FinalNodes = [Nodes for Nodes in processed_list if len(Nodes)>0]

print(f'Left with {len(FinalNodes)} graphs after performing a deep search for each graph.')
```

Anyway, you only need to run the notebook in order.

## Output
The output is a list of lists. Each one containing objects that says which graph it is (with the indices in the input), with which permutation of angle-types (with the index it has in the generated perms file), and some additional information used in the process. We added an aditional cell in the cases we obtained surviving graphs, to list and draw it.
