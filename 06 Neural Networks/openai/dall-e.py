#/usr/bin/env python
# OpenAI API demo

#######################################################################

from openai import OpenAI

# read API key from file
with open('PHYS395.key', 'r') as file:
    key = file.read().strip()

# establish client access
client = OpenAI(api_key=key)

# run query
response = client.images.generate(
    model="dall-e-3",
    prompt="a fluffy ragdoll cat",
    size="1792x1024",
    quality="hd",
    n=1,
)

# generated image URL
url = response.data[0].url; print(url)

#######################################################################

import requests
from PIL import Image
from io import BytesIO

# fetch content of URL
response = requests.get(url)
img = Image.open(BytesIO(response.content))

#######################################################################

import matplotlib.pyplot as plt

# render image
fig = plt.figure()
fig.gca().set_aspect('equal')
plt.imshow(img, interpolation='none')

plt.show()
