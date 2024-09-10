# Double-Pendulum
The code for the double pendulum sim

## Usage

### General

1. import it first 
```python 
import doublePendulum.py as dp
```
2. Run the `dp.sim` function.
3. To display the results you can use `dp.plott` or `dp.animate` functions.

### Sim

This function is the main function of the simulation, you should run it to generate the results of the simulation. It will reseve the initial parameters.  
example:  

```python
L1 = 3; L2 = 4; O1 = np.pi/2; O2 = 0
t, ang, cart, Eenergy = dp.sim(5, L1, L2, O1 , O2)
```

### Plott

This function creates static graphs of the results of the simulation  

- The function does not display the plot automatically. You should add `plt.show()`  

example:

```python
dp.plott(t,ang,cart,Eenergy)  
plt.show()
```

### Animation

Use this function to create and save an animation as mp4 file.  

> NOTE:
>> 1. This function takes only the time array and both postions of the pendulum.
>> 2. You must enter the size of the plot (recommended `l1 + l2`)
>> 3. If the name of the file given to the function already exists, it might be overwritten.
>> * For this function to work you must have _ffmpeg_ installed

#### **For more info check the functions's discription**  
  
---

## Requirements

### Necessary

- Python library _numpy_  
To install:

```powershell
pip install numpy
```

- Python library _scipy_  
To install:

```powershell
pip install scipy
```

### Optional

For `dp.plott` and `dp.animation` functions:

- Python library _matplotlib_  
To install:

```powershell
pip install matplotlib
```

For the `dp.animation` function you also need to install:

- Multimedia Framework _ffmpeg_. [Installation](https://ffmpeg.org/download.html)

**If you have difficulty installing the python libraries, and you have more then one user on the OS try add `-user`**
