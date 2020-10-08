#!/usr/bin/env python
# coding: utf-8

# # About this book
# 
# The github repository for this book is [https://github.com/ststewart/synestiabook2](https://github.com/ststewart/synestiabook2).<p>
#     
# This work is distributed under the MIT license:
#     
# ```
# MIT License
# 
# Copyright (c) 2020 Gigja Hollyday
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ```
# 
# 

# This book was compiled using the following python ecosystem:

# In[1]:


# Record the verions information for this book for posterity
import platform
print('python version: ',platform.python_version())
del platform
import matplotlib
print('matplotlib version: ', matplotlib.__version__)
del matplotlib
import numpy
print('numpy version: ', numpy.__version__)
del numpy
import scipy
print('scipy version: ', scipy.__version__)
del scipy
import rebound
print('rebound version: ', rebound.__version__)
del rebound
import ipywidgets
print('ipywidgets version: ', ipywidgets.__version__)
del ipywidgets
get_ipython().system('jupyter-book --version')


# In[ ]:




