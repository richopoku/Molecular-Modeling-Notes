class Atoms():

	def __init__(self):
		
# Calculates the distance between any two atoms with x, y, z coordinates and their corresponding charges


        def charge(self):
              ch = self.cha
              for item in ch:
                  for k,v in item.items():
                      charge_mul = v[0] * v[1] * self.charge_conv
                      self.charges.append({k:charge_mul})
       
        def distance(self):
               co = self.coo
               for item in co:
                   for k,v in item.items():
                       distan = np.sqrt(((v[0] - v[3])**2 + (v[1] - v[4])**2 + (v[2] - v[5])**2)) 
                       self.dist_po.append({k:distan})

#put all different dict of charges and distance into one dict......................................


        def dist_dict(self):
               dist_po = self.dist_po
               for element in dist_po:
                   for k,v in element.items():
                       self.newDict_dist_po[k] = v
                       self.newDict_dist_po1 = np.fromiter(self.newDict_dist_po.values(), dtype=float)
                       self.newDict_dist_po1 = np.array_split(self.newDict_dist_po1, 288)
                       

                       
        def charg_dict(self):
               All_charges_x = self.charges
               for element in All_charges_x:
                   for k,v in element.items():
                       self.newDict_All_charges[k] = v
                       self.newDict_All_charges1 = np.fromiter(self.newDict_All_charges.values(), dtype=float)
                       self.newDict_All_charges1 = np.split(self.newDict_All_charges1, [])
                       self.newDict_All_charges1 = np.asarray(self.newDict_All_charges1)
                       

#calculation of coulomb potential for all molecules       

        def coulomb_pot(self):
               self.coulombpot = self.newDict_All_charges1 / (self.newDict_dist_po1)
                    
               self.total_coul_potent = np.sum(self.coulombpot)
               self.coulomb_sum = np.sum(self.coulombpot, axis=1)
                


#combination of distance b/n 2 coordinates for atoms needed to compute the dispersion potential.....................

        def disper_coord(self):
             disper = self.disper
             for item in disper:
                  for k,v in item.items():
                      distance_disper = np.sqrt(((v[0] - v[3])**2 + (v[1] - v[4])**2 + (v[2] - v[5])**2)) 
                      self.dist_disper.append({k:distance_disper})

#put all different dict of distance into one dict and numpy array for the dispersion part.....

        def disper_indict(self):
             dist_disper = self.dist_disper
             for element in dist_disper:
                 for k,v in element.items():
                     self.newDict_dist_disper[k] = v
                     self.dist_dispersion = np.fromiter(self.newDict_dist_disper.values(), dtype=float)
                     self.dist_dispersion = np.array_split(self.dist_dispersion, 288)
                     self.dist_dispersion =  np.asarray(self.dist_dispersion)
                     
                     
#combination of distance b/n 2 coordinates for atoms needed to compute repulsion potential..................... 

        def repul_distext(self):
               repulse = self.repulsive_coord
               for item in repulse:
                   for k,v in item.items():
                       dista = np.sqrt(((v[0] - v[3])**2 + (v[1] - v[4])**2 + (v[2] - v[5])**2)) 
                       self.dist_repul.append({k:dista})

#put all different dict of distance into one dict and numpy array for the repulsion part.....

        def repul_indict(self):
             dist_repul = self.dist_repul
             for element in dist_repul:
                 for k,v in element.items():
                     self.newDict_dist_repul[k] = v
                     self.dist_repulsion = np.fromiter(self.newDict_dist_repul.values(), dtype=float)
                     self.dist_repulsion = np.array_split(self.dist_repulsion, 288)
                     self.dist_repulsion =  np.asarray(self.dist_repulsion)
                   
        
               
#calculation of dispersion and repulsion potential

        def damping_po(self):
              dist_dispersion = self.dist_dispersion
              self.damping = np.exp(-(1.28 * (self.R_m / dist_dispersion) - 1)**2)
              

         
        def dispersion_pot(self):
              damping = self.damping
              C = self.C_6
              dist_dispersion = self.dist_dispersion
              disper_pot = -(C / (dist_dispersion)**6 * np.exp(-(1.28 * (self.R_m / dist_dispersion) - 1)**2))         
              self.dispersionpotent = disper_pot 
              self.total_disp_potent = np.sum(self.dispersionpotent)
              self.dispersion_sum = np.sum(self.dispersionpotent, axis=1)
              


        def repulsion_pot(self):
              repul_pot = self.alpha * np.exp(-(1/0.529) * self.beta * (self.dist_repulsion))
              self.repulsionpotent = repul_pot
              self.total_repul_potent = np.sum(self.repulsionpotent)
              self.repulsion_sum = np.sum(self.repulsionpotent, axis=1)

        def total_disp_repul(self):
              self.total_disp_repul_pot = self.total_disp_potent + self.total_repul_potent


        #def error(self):
         #      self.water1 = ((self.theo1 - self.exp1)/ self.exp1) * 100
          #     self.waterall = ((self.theoall - self.expall)/ self.expall) * 100

        #def centermass1(self):
         #      center = self.newDict_dist_po1
          #     self.center_mass = np.sum(center, axis=1)
           #    self.center_mass = self.center_mass / 9

#Calculates the distance between two molecules based on their center of masses

        def centermass(self):
                center = self.Gwahcl
                dicGs = self.dicG
                self.center_mass = np.array_split(center, 288)
                self.center_mass  =  np.asarray(self.center_mass)
                self.newcenter_mass = dict(zip(dicGs,self.center_mass))
                
                for k,v in self.newcenter_mass.items():
                    dist_com = np.sqrt(((v[0] - v[3])**2 + (v[1] - v[4])**2 + (v[2] - v[5])**2)) 
                    self.dist_centcom.append({k:dist_com})
        def com(self):
               dist_cm = self.dist_centcom
               for element in dist_cm:
                   for k,v in element.items():
                       self.newDict_dist_cm[k] = v
                       self.newDict_dist_cm1 = np.fromiter(self.newDict_dist_cm.values(), dtype=float)
                       self.newDict_dist_cm1 = np.array_split(self.newDict_dist_cm1, 288)
                       self.newDict_dist_cm1 = np.asarray(self.newDict_dist_cm1)
                       self.newDict_dist_cm1 = np.concatenate(self.newDict_dist_cm1)
               
