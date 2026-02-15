import numpy as np
import spotipy
from spotipy.segmentation import get_masks

print("✅ Package imported successfully!")

# 1. Create a fake "Solar Disk" (10x10 pixels)
# The new code expects values between 0.0 and 2.0 (Standard Solar Intensity)
fake_sun = np.ones((10, 10))  # Background = 1.0 (Quiet Sun)
fake_sun[4:6, 4:6] = 0.2      # Spot = 0.2 (Dark Umbra)

print("\n--- Testing Segmentation ---")
print("Fake Sun Data (Center 4x4):\n", fake_sun[3:7, 3:7])

# 2. Test your function (Using the NEW arguments)
# We don't need to pass 'threshold' anymore.
# The function automatically normalizes 0.2 -> ~25 (which fits in the default umbra range 10-55)
masks = get_masks(fake_sun, cleanup=False)

# 3. Get the specific mask we want (Umbra)
umbra_mask = masks['umbra']

print("\nGenerated Umbra Mask (True = Spot detected):\n", umbra_mask[3:7, 3:7])

# 4. Verification
if umbra_mask[4, 4] == True and umbra_mask[0, 0] == False:
    print("\n✅ SUCCESS: The spot was detected correctly!")
else:
    print("\n❌ FAILURE: The spot was missed.")
