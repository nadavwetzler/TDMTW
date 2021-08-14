from obspy.clients.fdsn import Client, RoutingClient

def GetClient(fdsn):
    
    fdsn_names = fdsn.split(',')
    fdsn_l = len(fdsn_names)
    
    print('Using %d servers' % fdsn_l)
    servers = []

    for ii in range(fdsn_l):
        fdsn_name = fdsn_names[ii]
        if fdsn_name == 'GFZ':
            out = 0
            faultsT = 1
            client = Client(Client_def)
            CLIENT = Client_def
        elif fdsn_name == 'GFZ':
            CLIENT = "http://geofon.gfz-potsdam.de"
            out = 1
            faultsT = 2
            client = Client(CLIENT)
        elif fdsn_name == 'RC':
            client = RoutingClient("iris-federator")
            #client = RoutingClient("eida-routing")
            out = 1
            faultsT = 2
        elif fdsn_name == 'GNSr':
            CLIENT = "http://service-nrt.geonet.org.nz"
            # CLIENT = 'http://service.geonet.org.nz'
            out = 1
            faultsT = 2
            client = Client(CLIENT)
        elif fdsn_name == 'NOA':
            CLIENT = "http://eida.gein.noa.gr/"
            out = 1
            faultsT = 2
            client = Client(CLIENT)
        
        
        else:
            CLIENT = fdsn_name
            out = 1
            faultsT = 2
            client = Client(CLIENT)
        servers.append(client)
        print('Using FDSN: %s %s' % (client.base_url, CLIENT))
    client = servers

    return client, out, faultsT, fdsn_l
