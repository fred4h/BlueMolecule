#!/usr/bin/env python3
"""Write native/icons/BlueMolecule.ico (the Windows icon) from the rendered PNG icons in native/icons.

    python3 tools/make_ico.py        # needs Pillow; the .ico is kept in the repository, so builds do not need this

Every size is stored as a classic 32-bit BMP image with its AND mask (no PNG-compressed entries), the most
widely accepted form, including by the /win32icon option of the .NET Framework C# compiler that builds the
Windows launcher.
"""
import struct
from pathlib import Path

from PIL import Image

ICONS = Path(__file__).resolve().parent.parent / 'native' / 'icons'
SIZES = (16, 24, 32, 48, 64, 128, 256)


def rgba(size):
    p = ICONS / ('icon_%d.png' % size)
    src = p if p.exists() else ICONS / 'icon_256.png'
    im = Image.open(src).convert('RGBA')
    return im if im.size == (size, size) else im.resize((size, size), Image.LANCZOS)


def bmp_entry(im):
    s = im.size[0]
    px = im.load()
    xor = bytearray()
    for y in range(s - 1, -1, -1):                       # bottom-up rows, BGRA
        for x in range(s):
            r, g, b, a = px[x, y]
            xor += bytes((b, g, r, a))
    row = ((s + 31) // 32) * 4                           # 1 bit per pixel, rows padded to 32 bits
    mask = bytearray()
    for y in range(s - 1, -1, -1):
        bits = bytearray(row)
        for x in range(s):
            if px[x, y][3] == 0:                         # transparent pixel: set in the AND mask
                bits[x // 8] |= 0x80 >> (x % 8)
        mask += bits
    header = struct.pack('<IiiHHIIiiII', 40, s, 2 * s, 1, 32, 0, len(xor) + len(mask), 0, 0, 0, 0)
    return header + bytes(xor) + bytes(mask)


def main():
    images = [bmp_entry(rgba(s)) for s in SIZES]
    out = struct.pack('<HHH', 0, 1, len(images))
    offset = 6 + 16 * len(images)
    for s, data in zip(SIZES, images):
        out += struct.pack('<BBBBHHII', s % 256, s % 256, 0, 0, 1, 32, len(data), offset)
        offset += len(data)
    out += b''.join(images)
    (ICONS / 'BlueMolecule.ico').write_bytes(out)
    print('wrote %s (%d bytes, %d sizes)' % (ICONS / 'BlueMolecule.ico', len(out), len(SIZES)))


if __name__ == '__main__':
    main()
