from struct import unpack
import math
from PIL import Image

marker_mapping = {
    0xFFD8: "Start of Image",
    0xFFE0: "Application Default Header",
    0xFFDB: "Quantization Table",
    0xFFC0: "Start of Frame",
    0xFFC4: "Huffman Table",
    0xFFDA: "Start of Scan",
    0xFFD9: "End of Image",
}


def PrintMatrix(m):
    """
    A convenience function for printing matrices
    """
    for j in range(8):
        print("|", end="")
        for i in range(8):
            print("%d  |" % m[i + j * 8], end="\t")
        print()
    print()


def Clamp(col):
    """
    Makes sure col is between 0 and 255.
    """
    col = 255 if col > 255 else col
    col = 0 if col < 0 else col
    return int(col)


def ColorConversion(Y, Cr, Cb):
    """
    Converts Y, Cr and Cb to RGB color space
    """
    R = Cr * (2 - 2 * 0.299) + Y
    B = Cb * (2 - 2 * 0.114) + Y
    G = (Y - 0.114 * B - 0.299 * R) / 0.587
    return int(R + 128), int(G + 128), int(B + 128)
    # return (Clamp(R + 128), Clamp(G + 128), Clamp(B + 128))


def RemoveFF00(data):
    """
    Removes 0x00 after 0xff in the image scan section of JPEG
    """
    datapro = []
    i = 0
    while True:
        b, bnext = unpack("BB", data[i : i + 2])
        if b == 0xFF:
            if bnext != 0:
                break
            datapro.append(data[i])
            i += 2
        else:
            datapro.append(data[i])
            i += 1
    return datapro, i


def GetArray(type, l, length):
    """
    A convenience function for unpacking an array from bitstream
    """
    s = ""
    for i in range(length):
        s = s + type
    return list(unpack(s, l[:length]))


def DecodeNumber(code, bits):
    l = 2 ** (code - 1)
    if bits >= l:
        return bits
    else:
        return bits - (2 * l - 1)


class IDCT:
    """
    An inverse Discrete Cosine Transformation Class
    """

    def __init__(self):
        self.base = [0] * 64
        self.zigzag = [
            [0, 1, 5, 6, 14, 15, 27, 28],
            [2, 4, 7, 13, 16, 26, 29, 42],
            [3, 8, 12, 17, 25, 30, 41, 43],
            [9, 11, 18, 24, 31, 40, 44, 53],
            [10, 19, 23, 32, 39, 45, 52, 54],
            [20, 22, 33, 38, 46, 51, 55, 60],
            [21, 34, 37, 47, 50, 56, 59, 61],
            [35, 36, 48, 49, 57, 58, 62, 63],
        ]
        self.idct_precision = 8
        self.idct_table = [
            [
                (self.NormCoeff(u) * math.cos(((2.0 * x + 1.0) * u * math.pi) / 16.0))
                for x in range(self.idct_precision)
            ]
            for u in range(self.idct_precision)
        ]

    def NormCoeff(self, n):
        if n == 0:
            return 1.0 / math.sqrt(2.0)
        else:
            return 1.0

    def rearrange_using_zigzag(self):
        for x in range(8):
            for y in range(8):
                self.zigzag[x][y] = self.base[self.zigzag[x][y]]
        return self.zigzag

    def perform_IDCT(self):
        out = [list(range(8)) for i in range(8)]

        for x in range(8):
            for y in range(8):
                local_sum = 0
                for u in range(self.idct_precision):
                    for v in range(self.idct_precision):
                        local_sum += (
                            self.zigzag[v][u]
                            * self.idct_table[u][x]
                            * self.idct_table[v][y]
                        )
                out[y][x] = local_sum // 4
        self.base = out


class HuffmanTable:
    """
    A Huffman Table class
    """

    def __init__(self):
        self.root = []
        self.elements = []

    def BitsFromLengths(self, root, element, pos):
        if isinstance(root, list):
            if pos == 0:
                if len(root) < 2:
                    root.append(element)
                    return True
                return False
            for i in [0, 1]:
                if len(root) == i:
                    root.append([])
                if self.BitsFromLengths(root[i], element, pos - 1) == True:
                    return True
        return False

    def GetHuffmanBits(self, lengths, elements):
        self.elements = elements
        ii = 0
        for i in range(len(lengths)):
            for j in range(lengths[i]):
                self.BitsFromLengths(self.root, elements[ii], i)
                ii += 1

    def Find(self, st):
        r = self.root
        while isinstance(r, list):
            r = r[st.GetBit()]
        return r

    def GetCode(self, st):
        while True:
            res = self.Find(st)
            if res == 0:
                return 0
            elif res != -1:
                return res


class Stream:
    """
    A bit stream class with convenience methods
    """

    def __init__(self, data):
        self.data = data
        self.pos = 0

    def GetBit(self):
        b = self.data[self.pos >> 3]
        s = 7 - (self.pos & 0x7)
        self.pos += 1
        return (b >> s) & 1

    def GetBitN(self, l):
        val = 0
        for _ in range(l):
            val = val * 2 + self.GetBit()
        return val

    def len(self):
        return len(self.data)


class JPEG:
    """
    JPEG class for decoding a baseline encoded JPEG image
    """

    def __init__(self, image_file):
        # vertical and horizontal subsample factor for each component

        self.samp_factor = []
        # vertical and horizontal subsample interval for each component
        self.samp_interval = []
        # vertical and horizontal size for each component
        self.size = []

        # block dict for each component
        self.block = [{}, {}, {}]
        # image pixels
        self.pixel = {}

        self.huffman_tables = {}
        self.quant = {}
        self.quantMapping = []
        with open(image_file, "rb") as f:
            self.img_data = f.read()

    def DefineQuantizationTables(self, data):
        while len(data) > 0:
            (hdr,) = unpack("B", data[0:1])
            self.quant[hdr] = GetArray("B", data[1 : 1 + 64], 64)
            data = data[65:]

    def BuildMatrix(self, st, idx, quant, olddccoeff):
        i = IDCT()

        code = self.huffman_tables[0 + idx].GetCode(st)
        bits = st.GetBitN(code)
        dccoeff = DecodeNumber(code, bits) + olddccoeff

        i.base[0] = (dccoeff) * quant[0]
        l = 1
        while l < 64:
            code = self.huffman_tables[16 + idx].GetCode(st)
            if code == 0:
                break

            # The first part of the AC key_len
            # is the number of leading zeros
            if code > 15:
                l += code >> 4
                code = code & 0x0F

            bits = st.GetBitN(code)

            if l < 64:
                coeff = DecodeNumber(code, bits)
                i.base[l] = coeff * quant[l]
                l += 1

        i.rearrange_using_zigzag()
        i.perform_IDCT()

        return i.base, dccoeff

    def next(self, x, y, component):
        """
        :param x: horizontal coordinate of a block
        :param y: vertical coordinate of a block
        :param component: component idx
        :return: coordinate and component of the next block
        """
        if (x + 1)<self.size[component][1] and (x + 1)%self.samp_factor[component][1] != 0:
            # this block is not at the right edge of the image and mcu
            return x + 1, y, component
        if (y + 1)<self.size[component][0] and (y + 1)%self.samp_factor[component][0] != 0:
            # this block is not at the bottom of the image and mcu
            return x - x%self.samp_factor[component][1], y + 1, component
        if component != 2:
            # next component inside an mcu
            return (
                x // (self.samp_interval[component + 1][1] // self.samp_interval[component][1]),
                y // (self.samp_interval[component + 1][0] // self.samp_interval[component][0]),
                component + 1
            )
        if component == 2:
            if x + 1 == self.size[2][1] and y + 1 == self.size[2][0]:
                return None
            if x + 1 < self.size[2][1]:
                return (
                    (x + 1) * (self.samp_interval[2][1] // self.samp_interval[0][1]),
                    y * (self.samp_interval[2][0] // self.samp_interval[0][0]),
                    0
                )
            else:
                return (
                    0,
                    (y + 1) * (self.samp_interval[2][0] // self.samp_interval[0][0]),
                    0
                )

    def StartOfScan(self, data, hdrlen):
        data, lenchunk = RemoveFF00(data[hdrlen:])
        huffmanMapping = [0, 1, 1]

        st = Stream(data)

        # block descriptor: (x, y, component)
        block_descriptor = (0, 0, 0)
        prev_DC_coefficient = [0,0,0]
        
        while block_descriptor is not None:
            (
                self.block[block_descriptor[2]][(block_descriptor[0], block_descriptor[1])],
                prev_DC_coefficient[block_descriptor[2]]
            ) = self.BuildMatrix(
                st,
                huffmanMapping[block_descriptor[2]],
                self.quant[self.quantMapping[block_descriptor[2]]],
                prev_DC_coefficient[block_descriptor[2]]
            )
            block_descriptor = self.next(block_descriptor[0], block_descriptor[1], block_descriptor[2])

        return lenchunk + hdrlen

    def BaselineDCT(self, data):
        hdr, self.height, self.width, components = unpack(">BHHB", data[0:6])
        print("size %ix%i" % (self.width,  self.height))

        for i in range(components):
            id, samp, QtbId = unpack("BBB", data[6 + i * 3 : 9 + i * 3])
            self.quantMapping.append(QtbId)
            
            # get subsample info
            self.samp_factor.append([samp >> 4, samp & 0x0F])
        
        self.samp_interval = [
            [
                self.samp_factor[0][0] // self.samp_factor[i][0],
                self.samp_factor[0][1] // self.samp_factor[i][1]
            ]
            for i in range(len(self.samp_factor))
        ]

        # vertical and horizontal size in block for each component
        padded_height = math.ceil(self.height / (8 * self.samp_factor[0][0])) * (8 * self.samp_factor[0][0])
        padded_width = math.ceil(self.width / (8 * self.samp_factor[0][1])) * (8 * self.samp_factor[0][1])
        self.size = [
            [
                math.ceil(padded_height / (8 * self.samp_interval[i][0])),
                math.ceil(padded_width / (8 * self.samp_interval[i][1]))
            ]
            for i in range(len(self.samp_interval))
        ]

    def decodeHuffman(self, data):
        while len(data) > 0:
            offset = 0
            (header,) = unpack("B", data[offset : offset + 1])
            print(header, header & 0x0F, (header >> 4) & 0x0F)
            offset += 1

            lengths = GetArray("B", data[offset : offset + 16], 16)
            offset += 16

            elements = []
            for i in lengths:
                elements += GetArray("B", data[offset : offset + i], i)
                offset += i

            hf = HuffmanTable()
            hf.GetHuffmanBits(lengths, elements)
            self.huffman_tables[header] = hf
            data = data[offset:]

    def getPixel(self, x, y, component):
        block_coordinate = (
            x // (8 * self.samp_interval[component][1]),
            y // (8 * self.samp_interval[component][0])
        )
        return self.block[component][block_coordinate][
            y % (8 * self.samp_interval[component][0])//self.samp_interval[component][0]
        ][
            x % (8 * self.samp_interval[component][1])//self.samp_interval[component][1]
        ]

    def displayImage(self, image_name):
        with Image.new("RGB", (self.width, self.height), (0,0,0)) as im:
            for y in range(im.height):
                for x in range(im.width):
                    pixel_in_RGB = ColorConversion(
                            int(self.getPixel(x, y, 0)),
                            int(self.getPixel(x, y, 2)),
                            int(self.getPixel(x, y, 1))
                        )
                    im.putpixel(
                        (x, y),
                        (
                            pixel_in_RGB[0] & 0xFF,
                            pixel_in_RGB[1] & 0xFF,
                            pixel_in_RGB[2] & 0xFF
                        )
                    )
            im.show("{}.jpg".format(image_name))
            im.save("{}.bmp".format(image_name))

    def savePixels(self, file_out):
        # color space conversion
        for y in range(self.height):
            for x in range(self.width):
                self.pixel[(x, y)] = ColorConversion(
                            int(self.getPixel(x, y, 0)),
                            int(self.getPixel(x, y, 2)),
                            int(self.getPixel(x, y, 1))
                        )
        with open(file_out, "w") as f:
            for component in range(len(self.block)):
                f.write("component ID = {}\n".format(component))
                for y in range(self.height):
                    for x in range(self.width):
                        f.write("{:5}\t".format(self.pixel[(x, y)][component]))
                    f.write('\n')

    def decode(self):
        data = self.img_data
        while True:
            (marker,) = unpack(">H", data[0:2])
            print(marker_mapping.get(marker))
            if marker == 0xFFD8:
                data = data[2:]
            elif marker == 0xFFD9:
                return
            else:
                (len_chunk,) = unpack(">H", data[2:4])
                len_chunk += 2
                chunk = data[4:len_chunk]
                if marker == 0xFFC4:
                    self.decodeHuffman(chunk)
                elif marker == 0xFFDB:
                    self.DefineQuantizationTables(chunk)
                elif marker == 0xFFC0:
                    self.BaselineDCT(chunk)
                elif marker == 0xFFDA:
                    len_chunk = self.StartOfScan(data, len_chunk)
                data = data[len_chunk:]
            if len(data) == 0:
                break


if __name__ == "__main__":
    name = "img/ustc-banner.jpg"
    img = JPEG(name)
    img.decode()
    img.displayImage(name)
    #img.savePixels("{}.txt".format(name))