import pkg_resources
pkgs = ["numpy"]
installed_pkgs = {pkg.key for pkg in pkg_resources.working_set}
not_present = []
flag = 1
for p in pkgs:
    if p in installed_pkgs:
        flag=flag&1
    else:
        flag=flag&0
        not_present.append(p)
if not(flag):
    print "**** Some packages in python are not installed ****"
    print "\n".join([" * "+x for x in not_present])
    
