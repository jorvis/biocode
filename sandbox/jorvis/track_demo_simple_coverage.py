#!/usr/bin/env python3

"""
Writing a few demos for a new set of classes to export Biothings as graphical panels.

This is an example of a (very) simple CoverageTrack
"""

import math
import cairo
import random

track_width, track_height = 500, 50

surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, track_width, track_height)
ctx = cairo.Context (surface)


ctx.set_line_width(1)
ctx.set_source_rgb (1, 1, 1) # Solid color
ctx.paint()
ctx.set_source_rgb (0.5, 0.5, 0.5) # Solid color

# We'll use this to smoothen out the variation to make it look
#  more like a real coverage plot
last_height = 0

for i in range(0, track_width):
    ctx.move_to(i, track_height)
    height = random.randint(1 + int(last_height / 2), track_height - int(last_height / 2))
    ctx.line_to(i, track_height - height)
    last_height = height
    
ctx.stroke()
ctx.scale(track_width, track_height) # Normalizing the canvas
surface.write_to_png ("example.png") # Output to PNG
