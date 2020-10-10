htmlContent = "<!DOCTYPE html>\n\
<html>\n\t\
<head>\n\t\t\
<meta charset=\"UTF-8\">\n\t\t\
<script type=\"text/javascript\" src=\"data.js\"></script>\n\t\t\
<script type=\"text/javascript\" src=\"md.js\"></script>\n\t\
</head>\n\t\
<body>\n\t\t\
<h3>Module %(MODNO)s</h3>\n\t\t\
<div id=\"%(MODID)s\"></div>\n\t\t"

## Here goes the table
endContent = "<script>\n\t\t\
buildModule(%(MODNO)s);\n\t\t\
tableCreate(\"%(MODID)s\",%(MODNO)s);\n\t\t\
window.onbeforeunload = function(){\n\t\t\
window.scrollTo(0,0);\n\t\t\
}\n\t\t</script>\n\t\t\
<style>\n\t\t#%(MODID)s{\n\t\t\
position: relative;\n\t\t}\
</style>\n\t</body>\n</html>"

