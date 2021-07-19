from urllib.request import urlopen

def info_fetch(name):
    try:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + "%20".join(name.split(" ")) + "/property/CanonicalSMILES,IUPACName/CSV"
        ans = urlopen(url).read().decode("utf8")
        return ans
    except Exception as e:
        return e