import sys
from treenode import Node


def read_tree_string(instr):
    """
    this takes a newick string as instr
    and reads the string and makes the
    nodes and returns the root node
    """
    root = None
    index = 0
    nextchar = instr[index]
    start = True
    keepgoing = True
    curnode = None
    while keepgoing:
        if nextchar == "(":
            if start:
                root = Node()
                curnode = root
                start = False
            else:
                newnode = Node()
                curnode.add_child(newnode)
                curnode = newnode
        elif nextchar == "[":
            note = ""
            index += 1
            nextchar = instr[index]
            while True:
                if nextchar == "]":
                    break
                index += 1
                note += nextchar
                nextchar = instr[index]
            curnode.note = note
        elif nextchar == ',':
            curnode = curnode.parent
        elif nextchar == ")":
            curnode = curnode.parent
            index += 1
            nextchar = instr[index]
            name = ""
            while True:
                if (nextchar == ',' or nextchar == ')' or nextchar == ':'
                        or nextchar == ';' or nextchar == '['):
                    break
                name += nextchar
                if index < len(instr)-1:
                    index += 1
                    nextchar = instr[index]
                else:
                    keepgoing = False
                    break
            curnode.label = name
            index -= 1
        elif nextchar == ';':
            keepgoing = False
            break
        elif nextchar == ":":
            index += 1
            nextchar = instr[index]
            brlen = ""
            while True:
                if (nextchar == ',' or nextchar == ')' or nextchar == ':'
                        or nextchar == ';' or nextchar == '['):
                    break
                brlen += nextchar
                index += 1
                nextchar = instr[index]
            curnode.length = float(brlen)
            index -= 1
        elif nextchar == ' ':
            index += 1
            nextchar = instr[index]
        else:  # this is an external named node
            newnode = Node()
            curnode.add_child(newnode)
            curnode = newnode
            curnode.istip = True
            name = ""
            while True:
                if (nextchar == ',' or nextchar == ')' or nextchar == ':'
                        or nextchar == ';' or nextchar == '['):
                    break
                name += nextchar
                index += 1
                nextchar = instr[index]
            curnode.label = name
            index -= 1
        if index < len(instr) - 1:
            index += 1
        nextchar = instr[index]
    return root


def read_tree_file_iter(inf):
    info = open(inf, "r")
    for i in info:
        if len(i) > 2:
            yield read_tree_string(i.strip())
    info.close()


if __name__ == "__main__":
    for s in read_tree_file_iter(sys.argv[1]):
        s.number_tree()
        print(s.is_rooted())
        # print((s.get_newick_repr(True, True)+";"))
