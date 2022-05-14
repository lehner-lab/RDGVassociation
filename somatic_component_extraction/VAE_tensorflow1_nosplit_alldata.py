'''
Mischan Vali Pour 
Variational Autoencoder
adapted/modified from https://github.com/greenelab/tybalt/blob/master/tybalt_vae.ipynb
adapted/inspired from "Hand-On Maschine Learning with Scikit-Learn, Keras & Tensorflow" by Aurélien Géron December 2020
Outputs different performance parameters
Runs with tensorflow 1.15.5
'''

####import all important stuff
##import important modules
# Python ≥3.5 is required
import sys
import argparse #for parsing
print(sys.version) #print version

# Scikit-Learn ≥0.20 is required
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

##keras stuff
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Input, Dense, Lambda, Layer, Activation, BatchNormalization
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras import backend as K
from tensorflow.keras import metrics, optimizers, losses
from tensorflow.keras.callbacks import Callback

print(keras.__version__)
print(tf.__version__)


# Common imports
import numpy as np
import os

#pandas
import pandas as pd

# for pearson
from scipy.stats import pearsonr



####start getting info
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--learning_rate',
                    help='learning rate of the Adam optimizer')
parser.add_argument('-b', '--batch_size',
                    help='number of samples to include in each learning batch')
parser.add_argument('-e', '--epochs',
                    help='how many times to cycle through the full dataset')
parser.add_argument('-n', '--num_components',
                    help='latent space dimensionality (size)')
parser.add_argument('-t', '--dataset_training', 
                    help='training dataset, put in full name + direc')
parser.add_argument('-k', '--kappa', 
                    help='kappa, how strongly to linearly ramp up the KL loss after each epoch')
parser.add_argument('-d', '--depth', 
                    help='define whether there should be a layer between latent and input/ouput layer, if yes depth=2, else depth=1')
args = parser.parse_args()

# Set hyper parameters
learning_rate = float(args.learning_rate)
batch_size = int(args.batch_size)
epochs = int(args.epochs)
latent_dim = int(args.num_components)
dataset_training = args.dataset_training
kappa = float(args.kappa)
beta = K.variable(0) #KL loss weighting at first epoch

#decide on whether there should be a layer before input and latent layer
#make this layer
depth = int(args.depth)
hidden_dim =  latent_dim*2 #douple size of latent dimensions


###upload the training and test which has already been scaled and split
#convert to numpy array
train_df = pd.read_csv(dataset_training,sep='\t')
train_df_input = np.array(train_df.drop(['sample_short'], axis=1))
#print(train_df_input.shape)

# Set architecture dimensions
original_dim = train_df_input.shape[1]

# Random seed
seed = int(np.random.randint(low=0, high=10000, size=1))
np.random.seed(seed)

# Function for reparameterization trick to make model differentiable
# custom layer to sample the codings given mean and log_var
# samples from a normal distribution
def sampling(args):
    # Function with args required for Keras Lambda function
    z_mean, z_log_var = args
    # sample epsilon for a normal distribution with same shape as input columns aka phenotypes aka somatic features
    epsilon = K.random_normal(shape=tf.shape(z_mean), mean=0.,
                              stddev=1)
    # estimate latent codings by adding epsilon to mu and standard deviation which was learned
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z

## custom layer with loss
class CustomVariationalLayer(Layer):
    def __init__(self, **kwargs):
        self.is_placeholder = True
        super(CustomVariationalLayer, self).__init__(**kwargs)
    def vae_loss(self, x_input, x_decoded):
        reconstruction_loss = original_dim * \
                              metrics.mse(x_input, x_decoded) #using here mean squared error, multiplying with original dim since tensor already uses mean
        kl_loss = - 0.5 * K.sum(1 + latent_log_encoded -
                                K.square(latent_mean_encoded) -
                                K.exp(latent_log_encoded), axis=-1)
        return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss)) #combining reconstruction and KL loss and taking the mean
    def call(self, inputs):
        x = inputs[0]
        x_decoded = inputs[1]
        loss = self.vae_loss(x, x_decoded)
        self.add_loss(loss, inputs=inputs)
        return x
    

# implement warmup to slowly ramp up the KL loss (beta=0 is vanilla autoencoder and beta=1 full VAE)
# KL will be weighted by KL*beta
# modified code from https://github.com/keras-team/keras/issues/2595 from https://github.com/greenelab/tybalt/blob/master/tybalt_vae.ipynb
# idea of ladder VAE from https://arxiv.org/abs/1602.02282class WarmUpCallback(Callback):
class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa 
    def on_epoch_end(self, epoch, logs={}): #on_epoch_begin, # Behavior on each epoch
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)
    
    
#### ENCODER ####
#first dense layer to get mean and log var
#batch normalization
#relu as activation function
pheno_input = Input(shape=(original_dim, ))
z_shape = latent_dim

if depth == 1:
    latent_mean = Dense(latent_dim,
                        kernel_initializer='glorot_uniform')(pheno_input)
    latent_log_var = Dense(latent_dim,
                           kernel_initializer='glorot_uniform')(pheno_input)
elif depth == 2:
    hidden_dense = Dense(hidden_dim,
                         kernel_initializer='glorot_uniform')(pheno_input)
    hidden_dense_batchnorm = BatchNormalization()(hidden_dense)
    hidden_enc = Activation('relu')(hidden_dense_batchnorm)
    latent_mean = Dense(latent_dim,
                         kernel_initializer='glorot_uniform')(hidden_enc)
    latent_log_var = Dense(latent_dim,
                            kernel_initializer='glorot_uniform')(hidden_enc)


# batch normalization and activation for mu and log variance
latent_mean_batchnorm = BatchNormalization()(latent_mean)
latent_mean_encoded = Activation('relu')(latent_mean_batchnorm)

latent_log_batchnorm = BatchNormalization()(latent_log_var)
latent_log_encoded = Activation('relu')(latent_log_batchnorm)

# reparameterization
latent_codings = Lambda(sampling,
           output_shape=(z_shape, ))([latent_mean_encoded, latent_log_encoded])


#### DECODER ####
#need tanh as activation function since our output values are between -1 and 1
if depth == 1:
    decoder_to_reconstruct = Dense(original_dim,
                                   kernel_initializer='glorot_uniform',
                                   activation='tanh')
elif depth == 2:
    decoder_to_reconstruct = Sequential()
    decoder_to_reconstruct.add(Dense(hidden_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='relu',
                                     input_dim=latent_dim))
    decoder_to_reconstruct.add(Dense(original_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='tanh'))

phenos_reconstruct = decoder_to_reconstruct(latent_codings)

############### build up autoencoder ###############
adam = optimizers.Adam(lr=learning_rate)
vae_layer = CustomVariationalLayer()([pheno_input, phenos_reconstruct])
vae = Model(pheno_input, vae_layer)
vae.compile(optimizer=adam, loss=None, loss_weights=[beta])

############### fit Model ##########################
history = vae.fit(train_df_input,
                  shuffle=True,
                  epochs=epochs,
                  batch_size=batch_size,
                  validation_data=(train_df_input, None),
                  callbacks=[WarmUpCallback(beta, kappa)])

# evaluate final loss
training_loss= vae.evaluate(train_df_input)
#validation_loss= vae.evaluate(train_df_input)

print(training_loss)
#print(validation_loss)


############### check distribution of node activations #############
#encoder model
encoder = Model(pheno_input, latent_mean_encoded) #latent_log_encoded,latent_mean_encoded,latent_codings

#encoding on test set
encoded_train_df = encoder.predict_on_batch(train_df_input)
encoded_train_df = pd.DataFrame(encoded_train_df, columns= range(1, latent_dim+1))

sum_node_activity = encoded_train_df.sum(axis=0).sort_values(ascending=False)
#print(sum_node_activity)
sum_node_activity_mean = sum_node_activity.mean()


############### check Correlation with "golden set" independent components: UV, smoking, dHR, dMMR #############
##upload file
file_golden_IC = './input_files/TCGA_Hartwig_PCAWG_least45of56_allforGWAS_14832samples_smoking_UV_dHR_dMMR_ICs.txt'
golden_ICs = pd.read_csv(file_golden_IC,sep='\t')


##take the latent encoding matrix and this information to do pearson correlations
###### upload latent encodings of the mean
encoded_mean_train = encoder.predict_on_batch(train_df_input)
encoded_mean_train_df = pd.DataFrame(encoded_mean_train, columns= range(1, latent_dim+1))

#add sample short info
encoded_mean_dataset_df = encoded_mean_train_df.assign(sample_short= train_df.sample_short)


##################  output mean encoded layer  ##############################
output_mean_encoded_direc='./results/wholeData_meanEncoded'


output_latent_mean_encoded = os.path.join(output_mean_encoded_direc, "VAE_mse_plus_KL_" + str(latent_dim) + "_components_" + 
                               str(learning_rate) + "_lr_" + str(batch_size) + "_batchsize_" + 
                               str(epochs) + "_epochs_" + str(kappa) + "_kappa_" +
                               str(1) + "_beta_" + str(hidden_dim) + "_hiddenDim_" + str(depth) + "_depth_" + str(seed) + "_seed_latent_mean_encoded" + ".txt")
                                                                                      
encoded_mean_dataset_df.to_csv(output_latent_mean_encoded, sep='\t', index= False)

#############################################################################


##add golden ICs to encodings
encoded_encodings_IC_golden_set_ID = pd.merge(encoded_mean_dataset_df, golden_ICs, on='sample_short')

###### do pearson correlation between all columns except for sample_short ####
pearson_encodings_vs_ICs = encoded_encodings_IC_golden_set_ID.corr(method='pearson')

###estimate maximum correlation with the VAE_mean_encodings with each of the golden set components
UV_IC = max(pearson_encodings_vs_ICs.UV_IC.iloc[0:latent_dim])
Smoking_IC = max(pearson_encodings_vs_ICs.Smoking_IC.iloc[0:latent_dim])
dHR_IC = max(pearson_encodings_vs_ICs.dHR_IC.iloc[0:latent_dim])
dMMR_IC = max(pearson_encodings_vs_ICs.dMMR_IC.iloc[0:latent_dim])

##create df from it
max_correlation_golden_ICs = pd.DataFrame(np.array([[UV_IC,Smoking_IC,dHR_IC,dMMR_IC]]),
                                          columns=['UV_IC', 'Smoking_IC', 'dHR_IC', 'dMMR_IC'])
max_correlation_golden_ICs = np.absolute(max_correlation_golden_ICs) ##take absolute values/no minus

##get mean value
mean_pearson_golden_ICs = max_correlation_golden_ICs.mean(axis=1)
############# ############# ############# ############# ############# ############# ############# #############


############# pearson reconstruction vs input #################################################################
# encoding again
val_encoded = encoder.predict_on_batch(train_df_input)

####decoder generative model
decoder_input = Input(shape=(latent_dim, ))  # can generate from any sampled z vector
_x_decoded_mean = decoder_to_reconstruct(decoder_input)
decoder = Model(decoder_input, _x_decoded_mean)

# reconstruction
val_reconstructed = decoder.predict(val_encoded) 
val_reconstructed = pd.DataFrame(val_reconstructed, columns=train_df.drop(['sample_short'], axis=1).columns)
#print(val_reconstructed)

validation_df= pd.DataFrame(train_df_input, columns=train_df.drop(['sample_short'], axis=1).columns)

## check the mean pearson between input and reconstructed, pearson for each sample
r = [pearsonr(val_reconstructed.iloc[x, :],
              validation_df.iloc[x, :])[0] for x in range(val_reconstructed.shape[0])]
r_mean = np.mean(np.array(r))
#print(r_mean)


##################  output the results for a run ############################################################
output_parameter_model_run='./output/wholeData'

output = {'num_components': [latent_dim],
          'learning_rate': [learning_rate],
          'batch_size': [batch_size],
          'epochs': [epochs],
          'validation_loss': [0],
          'mean_sum_activity_mean_encoded_layer': [round(sum_node_activity_mean,4)],
          'r_mean_reconstruction': [round(r_mean,4)],
          'kappa': [kappa],
          'beta': [0],
          'UV_IC_max_cor': [round(UV_IC,4)],
          'Smoking_IC_max_cor': [round(Smoking_IC,4)],
          'dHR_IC_max_cor': [round(dHR_IC,4)],
          'dMMR_IC_max_cor': [round(dMMR_IC,4)],
          'mean_max_cor_golden_ICs': [round(mean_pearson_golden_ICs[0],4)],
          'depth': [depth],
          'hidden_layer_dim': [hidden_dim],
          'seed': [seed]
         }


output_df = pd.DataFrame(output, columns = ['num_components', 'learning_rate','batch_size', 'epochs',
                                            'validation_loss', 'mean_sum_activity_mean_encoded_layer',
                                            'r_mean_reconstruction','kappa','beta',
                                            'UV_IC_max_cor','Smoking_IC_max_cor','dHR_IC_max_cor',
                                            'dMMR_IC_max_cor','mean_max_cor_golden_ICs',
                                            'depth','hidden_layer_dim','seed'])

output_df_direc = os.path.join(output_parameter_model_run, "VAE_mse_plus_KL_" + str(latent_dim) + "_components_" + 
                               str(learning_rate) + "_lr_" + str(batch_size) + "_batchsize_" + 
                               str(epochs) + "_epochs_" + str(kappa) + "_kappa_" +
                               str(1) + "_beta_" + str(hidden_dim) + "_hiddenDim_" + str(depth) + "_depth_" + str(seed) + "_seed" + ".txt")
#save output
output_df.to_csv(output_df_direc, sep='\t', index= False)
##############################################################################################################

