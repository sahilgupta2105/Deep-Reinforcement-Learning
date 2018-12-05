from keras.layers import Input, Dense, Conv2D, UpSampling2D, AveragePooling2D
from keras.models import Model
from keras import backend as K
import numpy as np
from keras.callbacks import TensorBoard
from keras import regularizers
from keras.optimizers import Adam
from scipy.io import savemat

K.set_image_dim_ordering('tf')

def data_loader(x):
	n_steps = 500
	m = np.zeros((n_steps*x,900,2))
	m1 = np.zeros((n_steps*4,900,2))
	for i in range(1,x+1):
		for j in range(0,n_steps):
			temp_np = np.loadtxt('./train_data/'+str(i)+'/velocity_sampling/'+str(j)+'.txt',delimiter=',')
			m[j+(i-1)*n_steps,:,0] = temp_np[:,0]
			m[j+(i-1)*n_steps,:,1] = temp_np[:,1]

			if i<=4:
				temp_np = np.loadtxt('./test_data/'+str(i)+'/velocity_sampling/'+str(j)+'.txt',delimiter=',')
				m1[j+(i-1)*n_steps,:,0] = temp_np[:,0]
				m1[j+(i-1)*n_steps,:,1] = temp_np[:,1]

	# reshape the data
	x_train = np.zeros((n_steps*x,30,30,2))
	x_test = np.zeros((n_steps*4,30,30,2))

	for i in range(0,n_steps*x):
		x_train[i,:,:,0] = np.reshape(m[i,:,0],(30,30))
		x_train[i,:,:,1] = np.reshape(m[i,:,1],(30,30))

	for i in range(0,n_steps*4):
		x_test[i,:,:,0] = np.reshape(m1[i,:,0],(30,30))
		x_test[i,:,:,1] = np.reshape(m1[i,:,1],(30,30))
	savemat('ori_data.mat',{'o_m':x_train,'o_m1':x_test})

	return x_train, x_test

def build_autoencoder():
	
	init = 'RandomUniform' 
	act = 'relu'
	learn_rate = 0.001
	reg = regularizers.l2(0.00001)
	use_bias = False
	
	input_v = Input(shape=(30,30,2))
	conv1_1 = Conv2D(64, (3, 3), activation=act, padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(input_v)
	pool1 = AveragePooling2D((2,2),padding='same')(conv1_1)
	conv1_2 = Conv2D(64, (3,3),activation=act,padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(pool1)
	conv1_3 = Conv2D(64, (3,3),activation=act,padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(conv1_2)
	pool2 = AveragePooling2D((3,3),padding='same')(conv1_3)
	encoder = Dense(16,kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(pool2)
	decoder = Dense(16,kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(encoder)
	up1 = UpSampling2D((3,3))(decoder)
	conv2_1 = Conv2D(64,(3,3),activation=act,padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(up1)
	conv2_3 = Conv2D(64,(3,3),activation=act,padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(conv2_1)
	up2 = UpSampling2D((2,2))(conv2_3)
	conv2_2 = Conv2D(64,(3,3),activation=act,padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(up2)
	output_v = Conv2D(2,(3,3),activation='linear',padding='same',kernel_initializer=init,kernel_regularizer = reg,use_bias=use_bias)(conv2_2)

	autoencoder = Model(input=input_v,output=output_v)
	optimizer = Adam(lr=learn_rate)
	autoencoder.compile(optimizer='adam',loss='mse')

	return autoencoder;

x_train, x_test = data_loader(1)
print(x_train.shape)
print(x_test.shape)

model = build_autoencoder()
model.summary()
model.fit(x_train, x_train,
                epochs=50,
                batch_size=10,
                shuffle=True,
                validation_data = (x_test,x_test),
                callbacks=[TensorBoard(log_dir='/tmp/autoencoder12')])

# save model for using with TRPO
model.save('autoencoder.h5')

## writing data to check
decoded_train = model.predict(x_train)
decoded_test = model.predict(x_test)

savemat('train_pred.mat',{'m':decoded_train,'pred':decoded_test})