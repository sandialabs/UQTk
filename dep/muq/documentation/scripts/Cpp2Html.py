
import re
import markdown2

import sys


class TutorialDocument:

    def __init__(self, filename):

        # The regular expression that matches the special "/*** */" tutorial comment containing markdown
        self.exprString = '(\/\*\*\*.*?\*\/)'

        # Load the file
        with open(filename) as fin:
            self.fileContent = fin.read()

    def GetParts(self):
        """
        Splits the file contents into parts based on the tutorial comments.
        """
        pieces = re.split(self.exprString, self.fileContent, flags=re.S)
        return pieces

    def GetHeader(self):
        header ="""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css">
<link rel='stylesheet' href='https://fonts.googleapis.com/css?family=Titillium+Web:300' type='text/css'>
<link rel='stylesheet' href='https://fonts.googleapis.com/css?family=Arimo:400,700' type='text/css'>
<script src="https://code.jquery.com/jquery-2.1.4.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
<script src="https://cdn.rawgit.com/google/code-prettify/master/loader/run_prettify.js"></script>
<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {
    inlineMath: [['$','$'], ['\\(','\\)']],
    processEscapes: false
  }
});
</script>
<title>Example Code</title>
</head>

        <body>
        """

        return header

    def GetFooter(self):
        footer = """</body>"""
        return footer

    def FormatCode(self, string):
        if(string.isspace()):
            return '\n'
        else:
            code = self.StripLines(string.replace('<','&lt;').replace('>','&gt;'))
            output = '\n'
            if(len(code)>0):
                output = '\n<pre class="prettyprint lang-cpp">\n'
                output += code
                output += '\n</pre>\n'

            return output

    def FormatMarkdown(self, string):
        output = '\n'
        output += markdown2.markdown(string[4:-2].strip(), extras=["code-friendly"])
        output += '\n'

        return output

    def FormatCompleteCode(self, codePieces):
        stripPieces = []
        for piece in codePieces:
            if(len(piece)>0):
                stripPieces.append(piece.lstrip('\n').rstrip(' '))

        completeCode = self.StripLines(''.join(stripPieces))

        output = markdown2.markdown('#Complete Code')
        output += '\n<pre class="prettyprint lang-cpp">\n'
        output += completeCode
        output += '\n</pre>\n'

        return output


    def ToHTML(self):
        """
        Writes HTML containing a pretty version of the tutorial file
        """
        pieces = self.GetParts()

        output = self.GetHeader()

        codePieces = []

        # Loop through the pieces and format things for html
        for piece in pieces:
            isCode = True

            if(len(piece)>3):
                if(piece[0:4]=='/***'):
                    isCode = False

            if(isCode):
                output += self.FormatCode(piece)
                codePieces.append(piece)
            else:
                output += self.FormatMarkdown(piece)

        output += self.FormatCompleteCode(codePieces)

        output += self.GetFooter()

        return output

    def StripLines(self,s):
        lines = s.splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        while lines and not lines[-1].strip():
            lines.pop()
        return '\n'.join(lines)

if(len(sys.argv)!=3):
    print("\nERROR: Incorrect number of command line arguments to Cpp2Html.py")
    print("\nUSAGE:\n\tpython Cpp2Html.py <c++ inputfile> <html outputfile>\n\n")

assert(len(sys.argv)==3)
filename = sys.argv[1]

htmlString = TutorialDocument(filename).ToHTML()

with open(sys.argv[2],'w') as fout:
    fout.write(htmlString)
