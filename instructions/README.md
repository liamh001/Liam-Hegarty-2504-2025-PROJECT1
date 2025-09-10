# Generating project instructions

#### Weave.jl already installed
Assuming you have already installed `Weave.jl`.

1. Navigate to the current directory (`src/instructions`) in your terminal using `cd`.

2. Run from your terminal:
```
julia runweave.jl
```

3. Open the file `project1.html` in any browser - you can do so from the terminal via the command:
```
open project1.html
```

#### Installing Weave.jl

1. Open the Julia REPL in your terminal via:
```
julia
```

2. Change to package mode by entering `]`

3. Run the command:
```
pkg> add Weave
```
*Note: do **not** enter `pkg>` yourself - this should be part of the REPL interface after you press `]`.*