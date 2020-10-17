#!/usr/bin/env python
# coding: utf-8

# In[2]:


import requests

url = "https://github.com/ststewart/synestiabook2/tree/master/synestia-book/docs/syndef/"
fin = open('TE_Example03_Cool01_snapshot_10500_long', 'rb')
files = {'file': fin}
try:
    r = requests.post(url, files=files)


# In[ ]:





# In[7]:


import urllib.request
req = urllib.request.Request(url='https://github.com/ststewart/synestiabook2/tree/master/synestia-book/docs/syndef/TE_Example03_Cool01_snapshot_10500_long/')
f = urllib.request.urlopen(req)


# In[ ]:




