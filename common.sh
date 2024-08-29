# bash script to read the config file in yaml format
python3 -c 'import yaml, os, sys; sys.exit() if not os.path.exists("'$1'") else "ok"; f=open("'$1'"); config=yaml.safe_load(f); f.close(); print(config.get("'$2'",""))'
