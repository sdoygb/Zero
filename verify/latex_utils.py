#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility helpers for Markdown math in verification scripts."""
import re


BLOCK_MATH = re.compile(r"\$\$(.*?)\$\$", re.S)
INLINE_MATH = re.compile(r"\$([^$\n]+?)\$")


def canonical_math(text):
    """Normalize dollar math back to the backslash form used by early scripts."""
    text = BLOCK_MATH.sub(
        lambda match: "\\[\n" + match.group(1).strip("\n") + "\n\\]",
        text,
    )
    return INLINE_MATH.sub(
        lambda match: "$" + match.group(1) + "$",
        text,
    )
