### 1. Install the following MATLAB requirements and packages
- required MATLAB version: R2021B

### 2. Example
- Example code namely "main.m" is a MATLAB .m file that shows a demonstration of the frequency flattening method to reveal the hidden spirals. It takes a phase matrix as input and reveals the hidden spiral(s), if any, associated with it. Supporting functions are also provided.  

Input data expected in workspace:

- `DataCoordinate`
  - `DataCoordinate(1, :)` = trial number
  - `DataCoordinate(2, :)` = epoch start frame
  - `DataCoordinate(3, :)` = epoch end frame
- `Phases`
  - Phase time sequence
- `directory`
  - output folder

Pipeline order:

1. Extract 8x8 phase.
2. Interpolate to 29x29.
3. Find parameters for coupling function.
4. Check target stability.
5. Try continuation first.
6. If continuation is not accepted, try homotopy.

Parameters are not universal. You need to choose them yourself for the dataset and goal. Some epochs will still need hand tuning on numerical parameters.

Homotopy reversibility is not handled automatically here. For reversibility checks, use reverse homotopy, and for difficult epochs expect additional hand tuning.

`demo_two_cases.m`

Simple demo script for quick demonstration. Two saved example phases showing both methods and outcomes. 
  
