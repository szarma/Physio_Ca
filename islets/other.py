def create_download_command(glink, filename):
    from os import path
    fileid = glink.split("google.com")[1].strip("/").split("/")[2]
    command = r"""wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=FILEID" -O 'FILENAME' && rm -rf /tmp/cookies.txt"""
    command = command.replace("FILEID", fileid)
    command = command.replace("FILENAME", filename)
    if path.isfile(filename):
        print (f"{filename} already exist. please remove it and try again.")
        return ""
    return command
