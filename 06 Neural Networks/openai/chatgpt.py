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
completion = client.chat.completions.create(
    model="gpt-4o",
    messages=[
        {"role": "system", "content": "You are a helpful assistant."},
        {
            "role": "user",
            "content": "Write a haiku about recursion in programming."
        }
    ]
)

# generated message
msg = completion.choices[0].message; print(msg.content)
